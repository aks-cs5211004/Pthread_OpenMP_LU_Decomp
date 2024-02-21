#include "stdio.h"

__global__ void cuda_add(int a, int b, int *c){
	*c = a + b;
    printf("Hello from %d \n", threadIdx.x);
}

int main (){
	int a,b,c;
	int *dev_c;

	a=3;
	b=4;

	int driverVersion;
	int runtimeVersion;

	cudaDriverGetVersion(&driverVersion);
	cudaRuntimeGetVersion(&runtimeVersion);
	printf("%d, %d\n",driverVersion, runtimeVersion);


	cudaMalloc((void**)&dev_c, sizeof(int));
	cudaError_t err1 = cudaGetLastError();
	if(err1 != cudaSuccess)
		printf("Error %s\n",cudaGetErrorString(err1));
	
	cuda_add<<<1,256>>>(a,b,dev_c);
	cudaError_t err2 = cudaGetLastError();
	if(err2 != cudaSuccess)
		printf("Error %s\n",cudaGetErrorString(err2));
	
	cudaDeviceSynchronize();
	cudaError_t err3 = cudaGetLastError();
	if(err3 != cudaSuccess)
		printf("Error %s\n",cudaGetErrorString(err3));
	
	cudaMemcpy(&c, dev_c, sizeof(int), cudaMemcpyDeviceToHost);
	cudaError_t err4 = cudaGetLastError();
	if(err4 != cudaSuccess)
		printf("Error %s\n",cudaGetErrorString(err4));
	
	printf("%d + %d is %d\n",a,b,c);
	cudaFree(dev_c);
	return 0;
}
