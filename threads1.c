#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main(int argc, char* argv[]){
	double sum =0.0; double factor=0.0;
	int thread_count;
	thread_count=strtol(argv[1], NULL,10);
    int n;
	printf("Enter n");
	scanf("%d",&n);
#	pragma omp parallel for num_threads(thread_count) \
	reduction(+: sum) private(factor)
    for(int k=0;k<n;k++)
    {
		factor=(k%2==0)? 1.0:-1.0;
        sum+=factor/(2*k+1);
    }
	printf("%.40e\n",4*sum);
	return 0;
}




