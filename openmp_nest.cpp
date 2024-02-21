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

void luDecomposition(double *a, double* l, double* u, int* pi, int n, int total_threads)
{
	for(int k=0 ; k<n ; k++)
    {
        cout<<k<<endl;
        double max = 0.0;
        int kd = -1;
        for(int i=k; i<n ; i++)
        {
            if(max < abs(a[i*n + k]))
            {
                max = abs(a[i*n + k]);
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

        // omp_set_nested(1);
        #pragma omp parallel sections
        {
            #pragma omp section
            {
                double temp1;
                #pragma omp parallel for num_threads(total_threads/2) default(none) shared(k,kd,a,n) private(temp1)
                // swap a(k,:) and a(k',:)
                for(int i=0 ; i<n ; i++)
                {
                    temp1 = a[k*n + i];
                    a[k*n + i] = a[kd*n + i];
                    a[kd*n + i] = temp1;
                }
            }

            #pragma omp section
            {
                double temp2;
                #pragma omp parallel for num_threads(total_threads/2) default(none) shared(k,kd,l,n) private(temp2)
                // swap l(k,1:k-1) and l(k',1:k-1)
                for(int i=0 ; i<k ; i++)
                {
                    temp2 = l[k*n + i];
                    l[k*n + i] = l[kd*n + i];
                    l[kd*n + i] = temp2;
                }
            }
        }

        u[k*n + k] = a[k*n + k];

        for(int i=k+1 ; i<n ; i++)
        {
            l[i*n + k] = a[i*n + k] / u[k*n + k];
            u[k*n + i] = a[k*n + i];
        }

        // omp_set_nested(0);
        #pragma omp parallel for num_threads(total_threads) default(none) shared(k,l,u,a,n) collapse(2)
        for(int i=k ; i<n ; i++)
        {
            for(int j=k ; j<n ; j++)
            {
                a[i*n + j] = a[i*n + j] - l[i*n + k] * u[k*n + j];
            }
        }
    }

	return;
}

int main(int argc, char** argv)
{
	if (argc < 3)
	{
		cerr << "Usage: " << argv[0] << " <Size of the Matrix> <Number of Threads>" << endl;
		return 1;
	}

	int n = atoi(argv[1]);
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

	// // cout << "Completed Execution" << endl;

	// // cout << "Matrix a" << endl;
	// // print2dArray(a, n);
	// // cout << "Matrix l" << endl;
	// // print2dArray(l, n);
	// // cout << "Matrix u" << endl;
	// // print2dArray(u, n);

	// // cout << "=+=+=============" << endl;

    // auto cal_start = high_resolution_clock::now();
	// double *lu = (double*)malloc(sizeof(double)*(n*n));
	// matrixMult(l, u, lu, n);

	// // cout << "================" << endl;

	// // cout << "lu :" << endl;
	// // print2dArray(lu, n);

	// double* pa = (double*)malloc(sizeof(double)*(n*n));
	// double* p = (double*)malloc(sizeof(double)*(n*n));
	// p =  calPermutationMatrix(pi, n);
	// matrixMult(p, aCopy, pa, n);

	// // cout << "pa: " << endl;
	// // print2dArray(pa, n);

	// double* rm = calResidualMat(pa, lu, n);
	// cout << "L21 Norm: " << l2_1Norm(rm, n) << endl;
    // auto cal_end = high_resolution_clock::now();
    // duration<double> cal_time = cal_end - cal_start;
    // cout << "Calculation Time: " << cal_time.count() << " seconds." << endl;


	return 0;
}

void matrixMult(double* A, double* x, double* y, int n)
{
    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {   
            double temp = 0.0;
            #pragma omp parallel for reduction(+:temp) num_threads(omp_get_max_threads())
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

    #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
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

#pragma omp parallel for collapse(2) num_threads(omp_get_max_threads())
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
    #pragma omp parallel for reduction(+:norm) num_threads(omp_get_max_threads())
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