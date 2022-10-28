To run serial.c:
example:
gcc serial.c -lm
time ./a.out 7000 1


To run pthread.c:
example:
gcc pthread.c -lm -lpthread
time ./a.out 7000 16


To run openmp:
example:
gcc openmp.c -fopenmp -lm
export OMP_NUM_THREADS=16
time ./a.out 7000 16
