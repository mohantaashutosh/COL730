#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<pthread.h>
//variable definitions
int n, thread_count, effective_thread_count, k, size; // size of matrix n x n, number of threads 't'
long thread;
pthread_t* thread_handles;//initilizing pthread variables globally for fast access
pthread_mutex_t mutex;
double **P;
double **A;
double **A_original;
int *pi;
double **U;
double **L;
double **res_mat;
double max;
int k_dash;

//Call residual_matrix before calling this!!
double L21(){
    double ans=0;
    for(int j=0;j<n;j++){
        double col_sum=0;
        for(int i=0;i<n;i++){
            col_sum+=res_mat[i][j]*res_mat[i][j];
        }
        ans+=sqrt(abs(col_sum));
    }
    return ans;
}

//Call pi2P before running residual matrix!!
void residual_matrix(){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            res_mat[i][j]=0.0;
            for(int k=0;k<n;k++){
                res_mat[i][j]+=P[i][k]*A_original[k][j]-L[i][k]*U[k][j];
            }
        }
    }
}

void pi2P(){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            if(pi[i]==j){
                P[i][j]=1.0;
            }else{
                P[i][j]=0.0;
            }
        }
    }
}
double double_abs(double x){
    if(x<=0.0){
        return -x;
    }
    return x;
}

void print_mat(double **mat){
    printf("\n");
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            printf("%f ",mat[i][j]);
        }
        printf("\n");
    }
}
void print_LU(){
    print_mat(L);
    print_mat(U);
}
void verify(){
    pi2P();
    residual_matrix();
    printf("L21 Norm is: %f\n",L21());
}
void* update_LU(void* rank){
    long my_first_i = ((long)rank)*size + k + 1;
    long my_last_i = ((long)rank+1)*size + k;
    if(my_last_i >= n)
        my_last_i = n-1;
    while(my_first_i <= my_last_i){
        L[my_first_i][k] = A[my_first_i][k] * 1.0/U[k][k];
        U[k][my_first_i] = A[k][my_first_i];
        my_first_i++;
    }

}

void* update_A(void* rank){
    long my_first_i=((long)rank)*size + (long)k + 1;
    long my_last_i =((long)rank+1)*size + (long)k;
    if(my_last_i >= n)
        my_last_i = n-1;
    while(my_first_i <= my_last_i){
        for(int j = k+1; j < n; j++){
            A[my_first_i][j] -= L[my_first_i][k] * U[k][j];
        }
        my_first_i++;
    }
}
void* swap_k_kdash(void *rank){
    long my_first_i = ((long)rank)*size;
    long my_last_i = ((long)rank+1)*size-1;
    if(my_last_i>=k)
        my_last_i = k-1;
    double tmp_l;
    while(my_first_i <= my_last_i){
        tmp_l = L[k][my_first_i];
        L[k][my_first_i] = L[k_dash][my_first_i];
        L[k_dash][my_first_i] = tmp_l;
        my_first_i++;
    }
}
void* k_dash_find(void* rank){
    long my_first_i = ((long)rank)*size+k;
    long my_last_i = ((long)rank+1)*size+k-1;
    int my_k_dash = -1;
    double my_max = 0.0;
    if(my_last_i>=n)
        my_last_i = n-1;
    while(my_first_i<=my_last_i){
        if(max < double_abs(A[my_first_i][k])){
            my_max = double_abs(A[my_first_i][k]);
            my_k_dash = my_first_i;
        }
        my_first_i++;
    }
    pthread_mutex_lock(&mutex);
    if(my_max >= max){
        max = my_max;
        k_dash = my_k_dash;
    }
    pthread_mutex_unlock(&mutex);

}
void pthread(){
    thread_handles=malloc(32*sizeof(pthread_t));
    for(k = 0; k < n; k++){
        max = 0.0;
        effective_thread_count = thread_count;
        if(n-k < thread_count)
            effective_thread_count = n-k;
        size = (n-k)/effective_thread_count;
        if((n-k)%effective_thread_count!=0)
            size++;
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, k_dash_find, (void*) thread);
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);
        for(int i = k; i < n; i++){
            if(max < double_abs(A[i][k])){
                max = double_abs(A[i][k]);
                k_dash = i;
            }
        }
        if(max == 0){
            printf("singular matrix");
            return;
        }

        //swapping π[k] and π[k']
        int tmp_pi = pi[k];
        pi[k] = pi[k_dash];
        pi[k_dash] = tmp_pi;

        //swapping a(k,:) and a(k',:)
        double* tmp_a=A[k];
        A[k]=A[k_dash];
        A[k_dash]=tmp_a;

        //swapping l(k,1:k-1) and l(k',1:k-1)
        
        if(k!=0){
            effective_thread_count = thread_count;
            if(k < thread_count)
                effective_thread_count = k;
            size = k/effective_thread_count;
            if(k%effective_thread_count!=0)
                size++;
            for(thread = 0; thread < effective_thread_count; thread++)
                pthread_create(&thread_handles[thread], NULL, swap_k_kdash, (void*) thread);
            for(thread = 0; thread < effective_thread_count; thread++)
                pthread_join(thread_handles[thread], NULL);
        }


        U[k][k]=A[k][k];



        if(n-k-1<=0)
            continue;
        //updating LU
        effective_thread_count = thread_count;
        if(n-k-1<thread_count)
            effective_thread_count = n-k-1;
        size = (n-k-1)/effective_thread_count;
        if((n-k-1)%effective_thread_count!=0)
            size++;
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, update_LU, (void*) thread);
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);



        //updating A[i][j]
        effective_thread_count = thread_count;
        if(n-k-1<thread_count)
            effective_thread_count = n-k-1;
        size = (n-k-1)/effective_thread_count;
        if((n-k-1)%effective_thread_count!=0)
            size++;
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_create(&thread_handles[thread], NULL, update_A, (void*) thread);
        for(thread = 0; thread < effective_thread_count; thread++)
            pthread_join(thread_handles[thread], NULL);
    }
}


int main(int argc, char* argv[]){
    // To ensure the we get valid arguments
    pthread_mutex_init(&mutex, NULL);
    if(argc!=3){
        printf("invalid arguments");
        return -1;
    }

    n =  strtol(argv[1],NULL,10);
    thread_count = strtol(argv[2],NULL,10);
    if(thread_count > n)
        thread_count = n;



    P=malloc(n*sizeof(double*));
    A = malloc(n*sizeof(double*));
    A_original = malloc(n*sizeof(double*));
    pi = malloc(n*sizeof(int));
    U = malloc(n*sizeof(double*));
    L = malloc(n*sizeof(double*));
    res_mat = malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
        P[i]=malloc(n*sizeof(double));
        A[i]=malloc(n*sizeof(double));
        A_original[i]=malloc(n*sizeof(double));
        L[i]=malloc(n*sizeof(double));
        U[i]=malloc(n*sizeof(double));
        res_mat[i] = malloc(n*sizeof(double));
        pi[i] = i;
        for(int j=0;j<n;j++){
            A[i][j] = drand48();
            L[i][j] = 0.0;
            U[i][j] = 0.0;
            if(i==j){
                L[i][j] = 1.0;
            }
        }
    }

    // custom_test();
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            A_original[i][j] = A[i][j];
        }
    }


    clock_t t;
    t = clock();
    
    //LU Decomposition main code starts here..
    pthread();
    //LU Decomposition of main code over..
    free(thread_handles);
    pthread_mutex_destroy(&mutex);

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    // verify();
    printf("Time Take for Decomposition: %f\n",time_taken);
    return 0;
}