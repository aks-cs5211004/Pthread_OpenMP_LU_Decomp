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
vector<int> pi;
vector<vector<long double>> l, u, pim, residual; 
vector<vector<long double>> a, tempa;



void print(vector<vector<long double>>& vec){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            cout<<vec[i][j]<<" ";
        }
        cout<<endl;
    }
    cout<<endl;
}

void initil(){

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(1, 100);
    a.resize(N,vector<long double>(N,0));

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            a[i][j]=distribution(gen);
        }
    }
    pi.resize(N,0);
    tempa=a;
    l.resize(N,vector<long double>(N,0));
    pim.resize(N,vector<long double>(N,0));
    u.resize(N,vector<long double>(N,0)); 
    for(int i=0;i<N;i++) l[i][i]=1; 
    for(int i=0;i<N;i++) pi[i]=i;
}


struct arg_struct
{
    long thread;
    int kd;
} * args;



void* solvegreat(void* args){

        
        // PREREQS
        struct arg_struct* argsi=(struct arg_struct*)args;
        long my_rank=argsi->thread;
        cout<<"k:"<<k<<" kd:"<<kd<<" Thread:"<<my_rank<<endl;


        // // MAX RULES
        // long double local_max=INT64_MIN;
        
        // long noele=((N-k)*1ll)/thread_count;
        // long start=k+my_rank*noele;
        // long end=start+noele-1;
        // if(my_rank==thread_count-1){
        //     end+=((N-k)*1ll)%thread_count;
        // }
        // // MAX STARTS
        // for(int i0=start;i0<=end;i0++){
            
        // }




        // SWAP RULES
        long noele_0=((N)*1ll)/thread_count;
        long start_0=my_rank*noele_0;
        long end_0=start_0+noele_0-1;
        if(my_rank==thread_count-1){
            end_0+=((N)*1ll)%thread_count;
        }
        long noele_1=((k)*1ll)/thread_count;
        long start_1=my_rank*noele_1;
        long end_1=start_1+noele_1-1;
        if(my_rank==thread_count-1){
            end_1+=((k)*1ll)%thread_count;
        }
        // SWAP START
        cout<<start_0<<" "<<end_0<<endl;
        cout<<start_1<<" "<<end_1<<endl;
        for(int i1=start_0;i1<=end_0;i1++){
            swap(a[k][i1],a[kd][i1]);  
        }
        for(int i2=start_1;i2<=end_1;i2++){
            swap(l[k][i2],l[kd][i2]);  
        }
        // SWAP DONE
        // BARRIER_SWAP
        pthread_barrier_wait(&barrier1);
        // cout<<"barrier_1: "<<my_rank<<endl;


        // PART2 RULES
        u[k][k]=a[k][k];
        long noele_2=((N-k-1)*1ll)/thread_count;
        long start_2=k+1+my_rank*noele_2;
        long end_2=start_2+noele_2-1;
        if(my_rank==thread_count-1){
            end_2+=((N-k-1)*1ll)%thread_count;
        }
        // PART2 STARTS
        for(int i3=start_2;i3<=end_2;i3++){
            l[i3][k]=a[i3][k]/u[k][k];
            u[k][i3]=a[k][i3];
        }
        // PART2 ENDS
        // BARRIER_PART2
        pthread_barrier_wait(&barrier2);
        // cout<<"The Barrier 2: "<<my_rank<<endl;
        

        // PARTITION RULES
        long noele_3=((N-k-1)*1ll)/thread_count;
        long start_3=k+1+my_rank*noele_3;
        long end_3=start_3+noele_3-1;
        if(my_rank==thread_count-1){
            end_3+=((N-k-1)*1ll)%thread_count;
        }
        // PARTITION STARTS
        // CACHE HITS! ROW MAJOR!
        for(int i=start_3; i<=end_3;i++){
            double l_cache=l[i][k];
            for(int j=k+1;j<N;j++){
                a[i][j]=a[i][j]-l_cache*u[k][j];
            }
        }
        // PARTITION ENDS
        return 0;
}


void swapvec(vector<long double>& v1, vector<long double>& v2, int incl){
    for(int i=0;i<=incl;i++){
        swap(v1[i],v2[i]);
    }
}






void solve(){

    long thread;
    pthread_barrier_init(&barrier1, NULL, thread_count);
    pthread_barrier_init(&barrier2, NULL, thread_count);
    pthread_barrier_init(&barrier3, NULL, thread_count);
    thread_handles=(pthread_t*)malloc(thread_count*sizeof(pthread_t));
    pthread_mutex_init(&mutex_p, NULL);

    for(k=0;k<N;k++){
        long double max=0;
        kd=-1;

        for(int i=k;i<N;i++){
            if(max<abs(a[i][k])){
                max=abs(a[i][k]);
                kd=i;
            }
        }
        if(kd==-1) return;

        swap(pi[k],pi[kd]);

        cout<<"kd here"<<kd<<endl;

        // CREATE
        for(thread=0;thread<thread_count;thread++){
            args=(arg_struct*)malloc(sizeof(struct arg_struct)*1);
            args->thread=thread;
            args->kd=kd;
        pthread_create(&thread_handles[thread],NULL,solvegreat,(void*)args);
        }

        // JOIN
        for(thread=0; thread<thread_count;thread++){
        pthread_join(thread_handles[thread], NULL);
        }
    }
    free(thread_handles);
    pthread_mutex_destroy(&mutex_p);
    pthread_barrier_destroy(&barrier1);
    pthread_barrier_destroy(&barrier2);
    pthread_barrier_destroy(&barrier3);
}




void calcpim(){
    for(int i=0;i<N;i++){
        pim[i][pi[i]]=1.0;
    }

}

vector<vector<long double>> mul(vector<vector<long double>>& a, vector<vector<long double>>& b ){
    vector<vector<long double>> result(N,vector<long double>(N,0));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                result[i][j] +=a[i][k] * b[k][j];
            }
        }
    }
    return result;

}



void L2(){
    long double sumOfSquares =0.0;

    cout<<"The Great Mul"<<endl;
    vector<vector<long double>> PA=mul(pim,tempa);
    vector<vector<long double>> LU=mul(l,u);
    cout<<"The Great Mul Ends"<<endl;



    cout<<"AfterMath:------------"<<endl;
    print(pim);
    print(tempa);
    print(a);
    print(l);
    print(u);
    print(PA);
    print(LU);

    residual.resize(N,vector<long double>(N,0));
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            residual[i][j]=PA[i][j]-LU[i][j];
        }
    }
    for (const auto& row : residual) {
        for (auto element : row) {
            sumOfSquares += element * element;
        }
    }
    cout<<"The Norm Effect "<<sumOfSquares<<endl;
}

int main(int argc, char *argv[]){
    
    N=stoi(argv[1]);
    initil();
    cout<<"INIT DONE"<<endl;
    thread_count=strtol(argv[2], NULL, 10);
   
    auto start = high_resolution_clock::now();
    solve();
    auto time = duration_cast<seconds>(high_resolution_clock::now()- start);
    cout<<"Time : "<<time.count()<<endl;

    
    

    // INPUTz
    calcpim();
    L2();
    cout<<"Parallel Time : "<<time.count()<<endl;


}