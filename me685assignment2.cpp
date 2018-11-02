#include <bits/stdc++.h>
//#include "matplotlib-cpp/matplotlibcpp.h"

#define Tx0 300.00
#define Txl 100.00
#define C 0.00000005
using namespace std;

double relError(vector<double> x,vector<double> y){
    double diff=0;
    for(int i=0;i<x.size();i++){
        diff += (x[i]-y[i])*(x[i]-y[i]);
    }
    diff = sqrt(diff);
    cout<<"\n\t\trelative error: "<<diff;
    return diff;
}

vector< vector<double> > updateA(vector< vector<double> > x,vector<double> y,double del){
    for(int i=1;i<y.size()-2;i++){
        x[i][i]=-1*(2+4*C*del*del*y[i]*y[i]*y[i]);
    }
    return x;
}

vector<double> updateB(vector<double> x,vector<double> y,double del){
    for(int i=1;i<y.size()-2;i++){
        x[i]=-3*C*del*del*y[i]*y[i]*y[i]*y[i]+del*del;
    }
    return x;
}

vector<double> gaussSeidel(vector< vector<double> > A,vector<double> B,vector<double> T){
    vector<double> iter = T;
    int flag = 1,count=0;
    while(flag && count<50){
        count++;
        cout<<endl;
        for(int i=0;i<T.size();i++){
            T[i] = B[i] - inner_product(A[i].begin(),A[i].end(),T.begin(),0) + T[i]*A[i][i];
            T[i] = T[i]/A[i][i];
            cout<<T[i]<<"\t";
        }
        if(relError(T,iter)<0.001){
            flag=0;
            break;
        }
        iter = T;
    }
    return T;
}

vector<double> gaussSeidel2(vector< vector<double> > A,vector<double> B,vector<double> T){
    int flag=1,count=0;
    vector<double> x = T;
    while(flag && count<20){
        count++;
        for(int i=0;i<T.size();i++){
            T[i]=B[i]/A[i][i];
            for(int j=0;j<T.size();j++){
                if(j==i)
                    continue;
                T[i] = T[i] - ((A[i][j]/A[i][i])*x[j]);
                x[i] = T[i];
            }
            cout<<T[i]<<"\t";
        }
        cout<<"\n";
        if(relError(T,x)<0.001){
            flag=0;
            break;
        }
    }
    return T;
}


int main(int argc,char* argv[]) {
    int nodalPoints = atoi(argv[1]);
    vector<double> Tlast,B;
    vector<double> T(nodalPoints,0);
    vector< vector<double> > A(nodalPoints,vector<double>(nodalPoints,0));
    double delT = (Txl - Tx0)/(nodalPoints-1);
    double delX = 5.000/(nodalPoints-1);
    cout<<"delX:"<<delX<<endl;
    for(int i=0;i<nodalPoints;i++){
        Tlast.push_back(Tx0 + i*delT);
        B.push_back(-3*pow(delX,2)*C*pow(Tlast[i],4)+pow(delX,2));
        try{
            A.at(i).at(i-1) = 1;
            A.at(i).at(i+1) = 1;
            A.at(i).at(i) = -1*(2.00+4*pow(delX,2)*C*pow(Tlast[i],3));
        }
        catch(const std::out_of_range& oor){
            fill(A[i].begin(),A[i].end(),0);
            A[i][i]=1;
        }
    }
    B[0]=Tx0; B[nodalPoints-1]=Txl;
    int flag =1,count = 0;
    for(int i = 0; i<T.size();i++){
        cout<<T[i]<<"\t"<<Tlast[i]<<"\t"<<B[i]<<endl;
    }
    for(int i = 0; i<T.size();i++){
        cout<<endl;
        for(int j=0;j<T.size();j++)
            cout<<A[i][j]<<"\t";
    }
    cout<<"\nstarting gauss seidel iterations"<<endl;
    while(flag && count<20){
        count++;
        T = gaussSeidel(A,B,T);
        if(relError(T,Tlast)<0.01){
            flag=0;
            break;
        }
        Tlast = T;
        fill(T.begin(),T.end(),0);
        A = updateA(A,Tlast,delX);
        B = updateB(B,Tlast,delX);
    }
    for(int i=0;i<T.size();i++)
        cout<<Tlast[i]<<"\t";
    return 0;
}
