#include <bits/stdc++.h>
#include<iostream>
#include <set>
#include <algorithm> 
#include <vector>
#include <pthread.h>
// TEMPLATES AND FUNC->
typedef long long ll;
#define sz(a) ((int) (a).size())
using namespace std;
int N=0;
int t=0;
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
        cout<<"lol1"<<endl;
    pi.resize(N,0);
    tempa=a;
    l.resize(N,vector<long double>(N,0));
    pim.resize(N,vector<long double>(N,0));
    u.resize(N,vector<long double>(N,0)); 
    for(int i=0;i<N;i++) l[i][i]=1; 
    for(int i=0;i<N;i++) pi[i]=i;
}


void swapvec(vector<long double>& v1, vector<long double>& v2, int incl){
    for(int i=0;i<=incl;i++){
        swap(v1[i],v2[i]);
    }
}


void solve(){
    for(int k=0;k<N;k++){
        long double max=0;
        int kd=-1;

        // COLUMN major
        for(int i=k;i<N;i++){
            if(max<abs(a[i][k])){
                max=abs(a[i][k]);
                kd=i;
            }
        }
        if(kd==-1) return;
        swap(pi[k],pi[kd]);

        // row major
        swapvec(a[k],a[kd],N-1);
        swapvec(l[k],l[kd],k-1);
        u[k][k]=a[k][k];
        for(int i=k+1;i<N;i++){
            l[i][k]=a[i][k]/u[k][k];
            u[k][i]=a[k][i];
        }
        for(int i=k+1; i<N;i++){
            for(int j=k+1;j<N;j++){
                a[i][j]=a[i][j]-l[i][k]*u[k][j];
            }
        }
    }
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
    vector<vector<long double>> PA=mul(pim,tempa);
    vector<vector<long double>> LU=mul(l,u);
    print(pim);

    print(tempa);
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
    cout<<sumOfSquares<<endl;
}

int main(int argc, char *argv[]){
    
    // INPUTz
    N=stoi(argv[1]);
    initil();
    solve();
    calcpim();
    L2();



}