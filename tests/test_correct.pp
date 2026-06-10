#include<iostream>
using namespace std;
int main(int argc,char** argv)
{
	ifstream in1(argv[1]);
	ifstream in2(argv[2]);
	int m1,n1,m2,n2;
	in1>>m1>>n1;
	in2>>m2>>n2;
	if(m1!=m2||n1!=n2)
	{
		cout<<"Matrix dimensions dont match";
		return 0;
	}
	float *a=new float[m1*n1];
	float *b=new float[m2*n2];
	for(int i=0;i<m1;i++)
		for(int j=0;j<n1;j++)
			in1>>a[i*n1+j];
	for(int i=0;i<m2;i++)
		for(int j=0;j<n2;j++)
			in2>>b[i*n2+j];
	for(int i=0;i<m1;i++)
		for(int j=0;j<n1;j++)
			if(a[i*n1+j]!=b[i*n2+j])
			{
				cout<<"Matrices differ at position "<<i<<" "<<j<<endl;
				return 0;
			}
	cout<<"Matrices are equal";
}
