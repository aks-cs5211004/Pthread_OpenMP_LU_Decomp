#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
double Trap(double a, double b, int n);
double f(double al){
	return al*al;	
}
int main(int argc, char* argv[]){
	double global_result =0.0;
	double a,b;
	int n;
	int thread_count;
	thread_count=strtol(argv[1], NULL,10);
	printf("Enter a,b,n");
	scanf("%lf %lf %d", &a, &b, &n);
    double h=(b-a)/n;
    global_result=(f(a)+f(b))/2.0;
#	pragma omp parallel for num_threads(thread_count) \
	reduction(+: global_result)
    for(int i=1;i<=n-1;i++)
		global_result+=f(a+i*h);
global_result=h*global_result;
// 3.33337500000000e+08 -> 200 Parts 5 threads
	printf("%.14e\n", global_result);
	return 0;
}




