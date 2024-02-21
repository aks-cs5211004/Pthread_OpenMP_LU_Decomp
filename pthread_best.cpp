#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <time.h>
#include <cmath>
#include <omp.h>

using namespace std;
using namespace chrono;

void matrixMult(double* A, double* x, double* y, int n);
void print1DINTArray(int *arr, int size);
void print1dDOUBLEArray(double* arr, int size);
void print2dArray(double *arr, int size);
double* calResidualMat(double* A, double* B, int n);
double* calPermutationMatrix(int* pi, int n);
double l2_1Norm(double* matrix, int n);
double* copyArray(double* matrix, int n);
long thread_count;
pthread_t* thread_handles;
pthread_barrier_t barrier1;
double maxgl;
int kd, k, n;

void initialization(double *a, double *l, double *u, int* pi, int n)
{
	for (int i = 0; i < n; i++)
	{
		pi[i] = i;

		for (int j = 0; j < n; j++)
		{
			a[i*n + j] = ((double)(rand()%1000)) / 100.0;
			if (j > i)
			{
				u[i*n + j] = a[i*n + j];
				l[i*n + j] = 0.0;
			}
			else if (i == j)
			{
				u[i*n + j] = a[i*n + j];
				l[i*n + j] = 1.0;
			}
			else
			{
				u[i*n + j] = 0.0;
				l[i*n + j] = a[i*n + j];
			}
		}
	}
	return;
}


typedef struct thread_arguments {          
    double* a;
    double* l;
    double* u;
    int thread;
    int total;
} args;

void* solvegreat(void* arg){
        // ((args*)arg)->
        // PREREQS 
        int my_rank=((args*)arg)->thread;
        int thread_count=((args*)arg)->total;


        // PART2 RULES
        ((args*)arg)->u[k*n + k] = ((args*)arg)->a[k*n + k];
        long noele_2=((n-k-1)*1ll)/thread_count;
        long start_2=k+1+my_rank*noele_2;
        long end_2=start_2+noele_2-1;
        if(my_rank==thread_count-1){
            end_2+=((n-k-1)*1ll)%thread_count;
        }
        double ukk=((args*)arg)->u[k*n + k];
        for(int i=start_2 ; i<=end_2 ; i++)
        {
            ((args*)arg)->l[i*n + k] = ((args*)arg)->a[i*n + k] / ukk;
            ((args*)arg)->u[k*n + i] = ((args*)arg)->a[k*n + i];
        }
        pthread_barrier_wait(&barrier1);


        // PARTITION RULES        
        long noele_3=((n-k-1)*1ll)/thread_count;
        long start_3=k+1+my_rank*noele_3;
        long end_3=start_3+noele_3-1;        
        if(my_rank==thread_count-1){
            end_3+=((n-k-1)*1ll)%thread_count;
        }
        for(int i=start_3 ; i<=end_3 ; i++)
        {   
            double cache=((args*)arg)->l[i*n+k];
            for(int j=k+1 ; j<n ; j++)
            {
                ((args*)arg)->a[i*n+j] -= cache*((args*)arg)->u[k*n+j]; 
            }
        }

        return 0;
}



void luDecomposition(double *a, double* l, double* u, int* pi, int n, int thread_count)
{
    long thread;
    thread_handles=(pthread_t*)malloc(thread_count*sizeof(pthread_t));
    pthread_barrier_init(&barrier1, NULL, thread_count);
	for(k=0 ; k<n ; k++)
    {
        cout<<k<<endl;
        // print2dArray(a,n);
        maxgl = 0.0;
        kd = -1;
        for(int i=k; i<n ; i++)
        {
            if(maxgl < abs(a[i*n + k]))
            {
                maxgl = abs(a[i*n + k]);
                kd = i;
            }
        }
        
        if(kd == -1)
        {
            printf("\n\nSingular matrix ERROR\nProgram terminated with code 1\n\n");
            cout << k << '\n';
            return;
        }

        int temp0;
        temp0 = pi[k];
        pi[k] = pi[kd];
        pi[kd] = temp0;

        double temp1;
        for(int i=0 ; i<n ; i++)
        {
            temp1 = a[k*n + i];
            a[k*n + i] = a[kd*n + i];
            a[kd*n + i] = temp1;
        }
        double temp2;
        for(int i=0 ; i<k ; i++)
        {
            temp2 = l[k*n + i];
            l[k*n + i] = l[kd*n + i];
            l[kd*n + i] = temp2;
        }
        
        for(thread=0;thread<thread_count;thread++){
            args *in = (args *)malloc(sizeof(args));
	        in->a = a;							   
	        in->l = l;
	        in->u = u;
	        in->thread = thread;
            in->total = thread_count;
            pthread_create(&thread_handles[thread],NULL,solvegreat,(void*)(in));
        }
        // JOIN
        for(thread=0; thread<thread_count;thread++){
            pthread_join(thread_handles[thread], NULL);
        }
        
    }
    free(thread_handles);
    pthread_barrier_destroy(&barrier1);
    

	return;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cerr << "Usage: " << argv[0] << " <Size of the Matrix> <Number of Threads>" << endl;
		return 1;
	}

	n = atoi(argv[1]);
	int total_threads = atoi(argv[2]);

    cout << n << " " << total_threads << endl;

	double *a, *l, *u;
	a = (double*)malloc(sizeof(double)*(n*n));
	l = (double*)malloc(sizeof(double)*(n*n));
	u = (double*)malloc(sizeof(double)*(n*n));

	int *pi;
	pi = (int*)malloc(sizeof(int)*(n));

	srand(time(0));

	// cout << "Going Good" << endl;
	// Initialization
	auto start_init = high_resolution_clock::now();
	initialization(a, l, u, pi, n);
	auto end_init = high_resolution_clock::now();
	duration<double> elapsed_init = end_init - start_init;

	// cout << "Matrix a" << endl;
	// print2dArray(a, n);
	// cout << "Matrix l" << endl;
	// print2dArray(l, n);
	// cout << "Matrix u" << endl;
	// print2dArray(u, n);

	double* aCopy = copyArray(a, n);
	// cout << "Matrix aCopy" << endl;
	// print2dArray(aCopy, n);

	cout << "Initialization time: " << elapsed_init.count() << " seconds." << endl;

	// Sequential Algorithm
	auto start_seq = high_resolution_clock::now();
	luDecomposition(a, l, u, pi, n, total_threads);
	auto end_seq = high_resolution_clock::now();
	duration<double> elapsed_seq = end_seq - start_seq;

	cout << "Parallel Time: " << elapsed_seq.count() << " seconds." << endl;


    cout << "=+=+=============" << endl;
	double *lu = (double*)malloc(sizeof(double)*(n*n));
	matrixMult(l, u, lu, n);
	cout << "================" << endl;
	cout << "lu :" << endl;
	// print2dArray(lu, n);
	double* pa = (double*)malloc(sizeof(double)*(n*n));
	double* p = (double*)malloc(sizeof(double)*(n*n));
	p =  calPermutationMatrix(pi, n);
	matrixMult(p, aCopy, pa, n);
	// cout << "pa: " << endl;
	// print2dArray(pa, n);
	double* rm = calResidualMat(pa, lu, n);
	cout << "L21 Norm: " << l2_1Norm(rm, n) << endl;

	return 0;
}

void matrixMult(double* A, double* x, double* y, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            double temp = 0.0;
            for (int k = 0; k < n; k++)
            {
                temp += A[i*n + k] * x[k*n + j];
            }
            y[i*n + j] += temp;
        }
    }
}

void print1DINTArray(int *arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        cout << arr[i] << " ";
    } cout << endl;
    cout << endl;
    return;
}

void print1dDOUBLEArray(double* arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        cout << arr[i] << " ";
    } cout << endl;
    cout << endl;
    return;
}

void print2dArray(double *arr, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            cout << arr[i*size + j] << " ";
        } cout << endl;
    } cout << endl;

    return;
}

double* calResidualMat(double* A, double* B, int n)
{
    double* residual_matrix;
    residual_matrix = (double*)malloc(sizeof(double)*(n*n));
    for (int i = 0; i < n; i++)
    {   
        for (int j = 0; j < n; j++)
        {
            residual_matrix[i*n + j] = (A[i*n + j] - B[i*n + j]);
        }
    }
    return residual_matrix;
}

double* calPermutationMatrix(int* pi, int n)
{
    double* permutation_matrix  = (double*)malloc(sizeof(double)*(n*n));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (j == pi[i])
            {
                permutation_matrix[i*n + j] = 1.0;
            }
            else {
                permutation_matrix[i*n + j] = 0.0;
            }
        }
    }
    return permutation_matrix;
}

double l2_1Norm(double* matrix, int n)
{
    double norm = 0.0;

    // Compute L2 norm for each column and sum them up
    for (int j = 0; j < n; ++j) {
        double colNorm = 0.0;

        for (int i = 0; i < n; ++i) {
            colNorm += pow(matrix[i*n + j], 2);
        }

        norm += sqrt(colNorm);
    }

    return norm;
}

double* copyArray(double* matrix, int n)
{
	double *cpy = (double*)malloc(sizeof(double)*(n*n));

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cpy[i*n + j] = matrix[i*n + j];
		}
	}
	return cpy;
}