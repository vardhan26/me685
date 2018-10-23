#include<bits/stdc++.h>
#include "matplotlib-cpp/matplotlibcpp.h"
namespace plt = matplotlibcpp; 
namespace plt2 = matplotlibcpp; 
using namespace std; 


#define delT 1  //time stepsize
#define delX 0.01 //path stepsize
#define delN 0.001  //stepsize used during numerical integration
#define pi 3.14159265359
#define epsilon 0.00000 //for stopping criteria

//helper function for stopping criteria
double maxdiff(vector<double> a,vector<double> b){
	double max=0;
	for(int i=0;i<a.size();i++){
		if(fabs(a[i]-b[i])>max)
			max = fabs(a[i]-b[i]);
	}
	return max;
}

//tri-diagonal matrix algorithm
vector<double> TDMA(vector<double> a,vector<double> b,vector<double> c,vector<double> d){
	vector<double> x(d.size());
	c[0] = c[0]/b[0];
	d[0] = d[0]/b[0];
	for(int i=1;i<(int)d.size();i++){
		c[i] = c[i]/(b[i]-a[i]*c[i-1]);
		d[i] = (d[i]-d[i-1]*a[i])/(b[i]-c[i-1]*a[i]);
	}
	x[x.size()-1]=d[d.size()-1];
	for(int i=x.size()-2;i>=0;i--){
		x[i] = d[i]- c[i]*x[i+1];
	}
	return x;
}

//calculating "exact" solution using the analytical solution
vector<double> exact(float x,int maxtime){
	vector<double> T(maxtime,1-x);
	vector<double> C;
	C.push_back(0);
	vector<double> Tlast(maxtime);
	cout<<maxtime<<endl;
	int p=0;
	do{
		Tlast=T;
		for(int i=0;i<=1000;i++)
			C[C.size()-1] += 2*delN*((i*delN) - 1)*sin(((int)C.size())*pi*(i*delN));
		cout<<C[C.size()-1]<<"\n";
		for(int time=0;time<maxtime;time++){	
			T[time] += C[C.size()-1]*sin(C.size()*pi*x)*exp((-1)*time*pow(C.size()*pi,2));
			cout<<time<<"   "<<T[time]<<endl;}
		C.push_back(0);
		p++;
	}while(maxdiff(T,Tlast)>epsilon && p<100);
	for(int i=0;i<(int)T.size();i++){                                                                            
         cout<<T[i]<<" ";}
	return T;
}
// Driver program 
int main() 
{   int i;       
    double lambda = delT/(delX*delX);                                                                          
    vector<double> X(101);                                                                                    
    vector<float> path;
    for(int i =0; i<(int)X.size();i++){
	path.push_back(i*delX);	
	}                                                                                                          
    vector<double> a(101,-lambda);     // <- diagonals of                                                                       
    vector<double> b(101,1+2*lambda);  // <- Tri-Diagonal                                                                     
    vector<double> c(101,-lambda);     // <- Matrix                                                                     
    vector<double> d(101,0);           // <- Vector on RHS                                                                      
    d[0]=1;                            // values corrected at the boundaries                                                                       
    b[0]=1; b[100]=1;                   
    c[0]=0; c[100]=0;                                                                                         
    a[0]=0; a[100]=0;                                                                                         
    vector<vector<double>> T;                                                                                                                                                                 
    int j,count=0;
	vector<double> D = d;
	X = TDMA(a,b,c,d);     
	T.push_back(X);
	count++;                                                                                             
	while(maxdiff(D,X)>epsilon && count<=10000){
	D = X;
	d = D;
	                                                                                   
	X = TDMA(a,b,c,d);
 	T.push_back(X);
 	count++;		                                                                                                                                         
 	}
	cout<<T.size()<<" "<<count<<endl;
	
// plotting
	vector<double> E = exact(0.25,T.size());
	vector<int> Time;
	for(int i=0;i<T.size();i++)
		Time.push_back(i*delT);
	D.clear();
	for(int j=0;j<T.size();j++)
		D.push_back(T[j][25]);	
	plt::plot(Time,E,"r");
	plt::plot(Time,D,"r--");
	D.clear();	
	for(int j=0;j<T.size();j++)
		D.push_back(T[j][50]);
	vector<double> F = exact(0.50,T.size());	
	plt::plot(Time,F,"g");	
	plt::plot(Time,D,"g--");
	D.clear();	
	for(int j=0;j<T.size();j++)
		D.push_back(T[j][75]);	
	vector<double> G = exact(0.75,T.size());	
	plt::plot(Time,G,"b");
	plt::plot(Time,D,"b--");	
	plt::show();
	plt2::plot(path,X);
	plt2::show();
    return 0; 
} 


















