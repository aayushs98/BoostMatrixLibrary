#include<bits/stdc++.h>
#include"Matrix.hpp"

int n=3;

int main(){
	
	Matrix<double> m1(n,n), m2(n,n),m3(n,n);
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++) {
			m1(i,j) = i + j + 0.3;m2(i,j) = i + j + 0.3;
		}
	}

	m3 = m2;
	//m1 += m1 + m2 ; 

	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			cout << i << " " << j << " " << m1(i,j) << endl;
		}
	}

	cout<<endl;
	/*
	Matrix<int> m4(1,2);
	m4(0,0)=0;m4(0,1)=1;
	cout<<m4.Get_Column_Size()<<endl;

	cout<<endl;
*/
	//Matrix<int> m5(m4);
	//cout<<m5(0,1)<<endl;
	return 0;
}