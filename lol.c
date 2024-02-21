#include <stdio.h>
#include <stdlib.h>
#include <omp.h>


double Trap(double a, double b, int n);


int main(int argc, char* argv[]){
	double global_result =0.0;
	double a,b;
	int n;
	int thread_count;
	thread_count=strtol(argv[1], NULL,10);
	printf("Enter a,b,n");
	scanf("%lf %lf %d", &a, &b, &n);
#	pragma omp parallel num_threads(thread_count) \
	reduction(+: global_result)

{

		global_result+=Trap(a,b,n);
}

// 3.33337500000000e+08 -> 200 Parts 5 threads
			

	printf("%.14e\n", global_result);
	return 0;

}


double f(double al){
	return al*al;	
}

double Trap(double a, double b, int n){
	double h,x,my_result;
	double local_a, local_b;
	int i, local_n;


	int my_rank= omp_get_thread_num();
	int thread_count=omp_get_num_threads();
	h=(b-a)/n;
	local_n= n/thread_count;
	local_a= a+my_rank*local_n*h;
	local_b = local_a + local_n*h;
	my_result = (f(local_a) + f(local_b))/2.0;
	for (i = 1; i <= local_n - 1; i++) {
		x = local_a + i*h;
		my_result += f(x);
	}	
	my_result = my_result*h;
	printf("Hello from %d thread of %d which calculated %.14e \n", my_rank, thread_count, my_result);

	return my_result;

}





