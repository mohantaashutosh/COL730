#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
//variable definitions
int n, thread_count;
double **P;
double **A;
double **A_original;
int *pi;
double **U;
double **L;
double **res_mat;


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

void verify(){
    pi2P();
    residual_matrix();
    printf("L21 Norm is: %f\n",L21());
}


// void custom_test(){
//     A[0][0]=2.0;A[0][1]=4.0;A[0][2]=1.0;A[0][3]=-3;
//     A[1][0]=-1.0;A[1][1]=-2.0;A[1][2]=2.0;A[1][3]=4;
//     A[2][0]=4.0;A[2][1]=2.0;A[2][2]=-3.0;A[2][3]=5;
//     A[3][0]=5.0;A[3][1]=-4.0;A[3][2]=-3.0;A[3][3]=1;
// }

void serial(){
    for(int k = 0; k < n; k++){
        double max = 0.0;
        int k_dash;
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
        for(int j = 0; j < k; j++){
            double tmp_l = L[k][j];
            L[k][j] = L[k_dash][j];
            L[k_dash][j] = tmp_l;
        }


        U[k][k]=A[k][k];
        for(int i = k+1; i < n; i++){
            L[i][k] = A[i][k]*1.0/U[k][k];
            U[k][i] = A[k][i];
        }
        for(int i = k+1; i < n; i++){
            for(int j = k+1; j < n; j++){
                A[i][j] -= L[i][k]*U[k][j];
            }
        }
    }
}
int main(int argc, char* argv[]){
    // To ensure the we get valid arguments
    if(argc!=3){
        printf("invalid arguments");
        return -1;
    }

    n =  strtol(argv[1],NULL,10);
    thread_count = strtol(argv[2],NULL,10);



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
    serial();
    //LU Decomposition of main code over..

    t = clock() - t;
    double time_taken = ((double)t)/CLOCKS_PER_SEC;
    // verify();
    printf("Time Take for Decomposition: %f\n",time_taken);
    return 0;
}