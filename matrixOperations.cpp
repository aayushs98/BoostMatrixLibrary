#include<bits/stdc++.h>
#include"Matrix.hpp"

using namespace std;

int main()
{
	Matrix<float> m1(2,2);

	cout<<endl<<"m1 =="<<endl;

	m1(0,0) = 1.5; m1(0,1) = 2;
	m1(1,0) = 3.8; m1(1,1) = 4;
	
	cout<<m1(0,0)<<" "<<m1(0,1)<<endl;
	cout<<m1(1,0)<<" "<<m1(1,1)<<endl;

	cout<<endl<<"m1 += m1 =="<<endl;
	
	m1+=m1;
	
	cout<<m1(0,0)<<" "<<m1(0,1)<<endl;
	cout<<m1(1,0)<<" "<<m1(1,1)<<endl;


	cout<<endl<<"m2 =="<<endl;

	Matrix<float> m2(2,2);

	m2(0,0) = 1.5; m2(0,1) = 2.5;
	m2(1,0) = 3.5; m2(1,1) = 4.5;
	
	cout<<m2(0,0)<<" "<<m2(0,1)<<endl;
	cout<<m2(1,0)<<" "<<m2(1,1)<<endl;

	Matrix<float> m3(2,2);

	cout<<endl<<"m1 + m2 =="<<endl;

	m3=m2+m1;

	cout<<m3(0,0)<<" "<<m3(0,1)<<endl;
	cout<<m3(1,0)<<" "<<m3(1,1)<<endl;

	cout<<endl<<"m2 * m3 =="<<endl;

	Matrix<float> m4(1,2);

	m4(0,0)=5.5;
	m4(0,1)=10.75;
	
	cout<<endl<<"m4 =="<<endl;

	cout<<m4(0,0)<<" "<<m4(0,1)<<endl;
	
	Matrix<float> m5(1,2);
	
	cout<<endl<<"m4 * m2 =="<<endl;

	m5=m4*m2;

	cout<<m5(0,0)<<" "<<m5(0,1)<<endl;
	
	cout<<endl<<"m2 *= m2 =="<<endl;

	m2*=m2;

	cout<<m2(0,0)<<" "<<m2(0,1)<<endl;
	cout<<m2(1,0)<<" "<<m2(1,1)<<endl;

	return 0;
}
