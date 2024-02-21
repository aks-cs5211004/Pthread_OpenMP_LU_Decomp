#include <bits/stdc++.h>
#include<iostream>
#include <set>
#include <algorithm> 
#include <vector>
#include <pthread.h>
#include <chrono>
using namespace std::chrono;
// TEMPLATES AND FUNC->
typedef long long ll;
#define sz(a) ((int) (a).size())
using namespace std;
int N=0;
long thread_count;
int k=0,kd=0; 
pthread_t* thread_handles;
pthread_barrier_t barrier1, barrier2, barrier3;
pthread_mutex_t mutex_p;
long double maxgl=0;
const unsigned long long len=25000000;
// 64000000
int* pi;
long double* l; 
long double* u; 
long double* tempa; 
long double* a;




void print(long double* vec){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            cout<<vec[i+j*N]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void initil(){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(1, 100);
    l = (long double*)malloc((N*N) * sizeof(long double));
    u =(long double*)malloc((N*N) * sizeof(long double));
    a = (long double*)malloc((N*N) * sizeof(long double));
    tempa = (long double*)malloc((N*N) * sizeof(long double));
    pi = (int*)malloc((N) * sizeof(int));

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            a[i+j*N]=distribution(gen);
            tempa[i+j*N]=a[i+j*N];
        }
    }
    for(int i=0;i<N;i++) l[i+N*i]=1; 
    for(int i=0;i<N;i++) pi[i]=i;
}




void* solvegreat(void* thread){

        
        // PREREQS
        long my_rank=(long)thread;
        


        
        // // MAX RULES
        // long double local_max=0;
        // int local_kd=-1;
        // long noele=((N-k)*1ll)/thread_count;
        // long start=k+my_rank*noele;
        // long end=start+noele-1;
        // if(my_rank==thread_count-1){
        //     end+=((N-k)*1ll)%thread_count;
        // }
        // for(int i0=start;i0<=end;i0++){
        //     if(local_max<abs(a[i0+N*k])){
        //         local_max=abs(a[i0+N*k]);
        //         local_kd=i0;
        //     }   
        // }
        // if(local_max>maxgl){
        //     pthread_mutex_lock(&mutex_p);
        //     maxgl=local_max;
        //     kd=local_kd;
        //     pthread_mutex_unlock(&mutex_p);
        // }
        // pthread_barrier_wait(&barrier1);
        // if(kd==-1) return 0;


        // // SWAP RULES
        // long noele_0=((N)*1ll)/thread_count;
        // long start_0=my_rank*noele_0;
        // long end_0=start_0+noele_0-1;
        // if(my_rank==thread_count-1){
        //     end_0+=((N)*1ll)%thread_count;
        // }
        // long noele_1=((k)*1ll)/thread_count;
        // long start_1=my_rank*noele_1;
        // long end_1=start_1+noele_1-1;
        // if(my_rank==thread_count-1){
        //     end_1+=((k)*1ll)%thread_count;
        // }
        // for(int i1=start_0;i1<=end_0;i1++){
        //     swap(a[k+i1*N],a[kd+N*i1]);  
        // }
        // for(int i2=start_1;i2<=end_1;i2++){
        //     swap(l[k+N*i2],l[kd+N*i2]);  
        // }
        // pthread_barrier_wait(&barrier2);
        if(my_rank==0)
            swap(pi[k],pi[kd]);


        // PART2 RULES
        
        u[k+k*N]=a[k+N*k];
        long noele_2=((N-k-1)*1ll)/thread_count;
        long start_2=k+1+my_rank*noele_2;
        long end_2=start_2+noele_2-1;
        if(my_rank==thread_count-1){
            end_2+=((N-k-1)*1ll)%thread_count;
        }
        for(int i3=start_2;i3<=end_2;i3++){
            l[i3+N*k]=a[i3+N*k]/u[k+N*k];
            u[k+N*i3]=a[k+N*i3];
        }
        pthread_barrier_wait(&barrier3);
        
        

        // PARTITION RULES
        long noele_3=((N-k-1)*1ll)/thread_count;
        long start_3=k+1+my_rank*noele_3;
        long end_3=start_3+noele_3-1;
        if(my_rank==thread_count-1){
            end_3+=((N-k-1)*1ll)%thread_count;
        }

        // cout<<my_rank<<"chunk size"<<(start-k-1-(end-k-1))*(start-k-1-(end-k-1))<<endl;
        for(int i=start_3; i<=end_3;i++){
            long double l_cache=l[i+N*k];
            for(int j=k+1;j<N;j++){
                a[i+N*j]-=l_cache*u[k+N*j];
            }
        }

        return 0;
}



void solve(){

    long thread;
    pthread_barrier_init(&barrier1, NULL, thread_count);
    pthread_barrier_init(&barrier2, NULL, thread_count);
    pthread_barrier_init(&barrier3, NULL, thread_count);
    
    thread_handles=(pthread_t*)malloc(thread_count*sizeof(pthread_t));
    pthread_mutex_init(&mutex_p, NULL);
    
    auto start = high_resolution_clock::now();
    for(k=0;k<N;k++){
        cout<<"outer_loop_var="<<k<<endl;
        maxgl=0;
        kd=-1;

        // SEQ PART
        for(int i=k;i<N;i++){
            if(maxgl<abs(a[i+N*k])){
                maxgl=abs(a[i+N*k]);
                kd=i;
            }
        }
        if(kd==-1) return;
        swap(pi[k],pi[kd]);
        for(int i=0;i<=N-1;i++){
            swap(a[k+N*i],a[kd+N*i]);
        }
        for(int i=0;i<=k-1;i++){
            swap(l[k+N*i],l[kd+N*i]);
        }
        // u[k+N*k]=a[k+N*k];
        // for(int i=k+1;i<N;i++){
        //     l[i+N*k]=a[i+N*k]/u[k+N*k];
        //     u[k+N*i]=a[k+N*i];
        // }
        // SEQ PART END




        // CREATE
        for(thread=0;thread<thread_count;thread++){
        pthread_create(&thread_handles[thread],NULL,solvegreat,(void*)thread);
        }
        // JOIN
        for(thread=0; thread<thread_count;thread++){
        pthread_join(thread_handles[thread], NULL);
        }
    }
    auto time = duration_cast<milliseconds>(high_resolution_clock::now()- start);
    cout<<"Time for Parallelised Code (in s) : "<<(long double)time.count()/1000<<endl;
    free(thread_handles);
    pthread_mutex_destroy(&mutex_p);
    pthread_barrier_destroy(&barrier1);
    pthread_barrier_destroy(&barrier2);
    pthread_barrier_destroy(&barrier3);
}




long double* calcpim(){
    long double* pim = new long double[N*N];
    for(int i=0;i<N;i++){
        pim[i+N*pi[i]]=1.0;
    }
    return pim;

}

long double * mul(long double*  a, long double * b ){
    long double* result = new long double[N*N];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                result[i+j*N] +=a[i+N*k] * b[k+N*j];
            }
        }
    }
    return result;

}



void L2(){
    long double sumOfSquares =0.0;

    cout<<"Calculating PA-LU"<<endl;
    long double * pim = calcpim();
    long double *  PA=mul(pim,tempa);
    long double* LU=mul(l,u);

    cout<<"Results:----------"<<endl;
    // cout<<"PI="<<endl;
    // print(pim);
    // cout<<"A="<<endl;
    // print(tempa);
    // cout<<"A="<<endl;
    // print(a);
    // cout<<"L="<<endl;
    // print(l);
    // cout<<"U="<<endl;
    // print(u);
    // cout<<"PA="<<endl;
    // print(PA);
    // cout<<"LU="<<endl;
    // print(LU);

    long double residual[N*N];
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            residual[i+N*j]=PA[i+N*j]-LU[i+N*j];
        }
    }
    for (int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            sumOfSquares += residual[i+N*j]*residual[i+N*j];
        }
    }
    cout<<"L2,1 Norm:  "<<sumOfSquares<<endl;
}

int main(int argc, char *argv[]){
    
    N=stoi(argv[1]);
    initil();
    cout<<"initialising....."<<endl;
    thread_count=strtol(argv[2], NULL, 10);
    solve();
    L2();


}