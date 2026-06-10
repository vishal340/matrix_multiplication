#include<iostream>
#include<fstream>
#include<cmath>
#include<cstring>
#include<sys/time.h>

int threshold=32;

using namespace std;


inline void add(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2)
{
for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++)
    C[i*jump2+j]=A[i*jump+j]+B[i*jump1+j];
}

inline void atomic_add(float* A,int jump,float* B,int jump1,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]+=A[i*jump+j];
}

inline void subtract(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      C[i*jump2+j]=A[i*jump+j]-B[i*jump1+j];
}

inline void atomic_subtract(float* A,int jump,float* B,int jump1,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]-=A[i*jump+j];
}

void multiply(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2,int n3)
{
for(int j=0;j<n1;j+=2)
  for(int i=0;i<n2;i+=2)
    for(int k=0;k<n3;k++)
    {
      C[j*jump2+k]+=A[i+j*jump]*B[i*jump1+k];
      C[j*jump2+k]+=A[i+1+j*jump]*B[(i+1)*jump1+k];
      C[(j+1)*jump2+k]+=A[i+(j+1)*jump]*B[i*jump1+k];
      C[(j+1)*jump2+k]+=A[(i+1)+(j+1)*jump]*B[(i+1)*jump1+k];
    }
}

inline void nonzero_C(float* C,int jump2,int m1,int m3)
{
  atomic_add(C+m3/2,jump2,C,jump2,m1/2,m3/2);
  atomic_add(C+jump2*m1/2,jump2,C+jump2*m1/2+m3/2,jump2,m1/2,m3/2);
  atomic_subtract(C,jump2,C+jump2*m1/2+m3/2,jump2,m1/2,m3/2);
}

void strassen(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2,int n3,int iter,bool flag)
{
  if(iter==1)
    multiply(A,jump,B,jump1,C,jump2,n1,n2,n3);
  else
  {
    iter/=2;
    int m1=n1/2,m2=n2/2,m3=n3/2;
    float* temp1=new float[m1*m2];
    float* temp2=new float[m2*m3];

    //M1
    add(A,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);   //temp1=A11+A22
    add(B,jump1,B+jump1*m2+m3,jump1,temp2,m3,m2,m3);   //temp2=B11+B22
    if(iter>1 && !flag)
      nonzero_C(C,jump2,m1,m3);
    strassen(temp1,m2,temp2,m3,C,jump2,m1,m2,m3,iter,flag);   //C11+=temp1*temp2
    if(flag)
    {
      for(int i=0;i<m1;i++)
        for(int j=0;j<m3;j++)
          C[jump2*m1+m3+i*jump2+j]=C[i*jump2+j];  //C22+=C11
    }
    else
      atomic_add(C,jump2,C+jump2*m1+m3,jump2,m1,m3);  //C22+=C11

    //M2
    add(A+jump*m1,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);  //temp1=A21+A22
    if(iter>1 && !flag)
      nonzero_C(C+jump2*m1,jump2,m1,m3);
    strassen(temp1,m2,B,jump1,C+jump2*m1,jump2,m1,m2,m3,iter,flag);  //C21+=temp1*B11

    //M5
    add(A,jump,A+m2,jump,temp1,m2,m1,m2);  //temp1=A11+A12
    if(iter>1 && !flag)
      nonzero_C(C+m3,jump2,m1,m3);
    strassen(temp1,m2,B+jump1*m2+m3,jump1,C+m3,jump2,m1,m2,m3,iter,flag);  //C12+=temp1*B22
    atomic_subtract(C+m3,jump2,C,jump2,m1,m3);  //C11-=C12

    //M6
    subtract(A +jump*m1,jump,A ,jump,temp1,m2,m1,m2);   //temp1=A21-A11
    add(B,jump1,B+m3,jump1,temp2,m3,m2,m3);    //temp2=B11+B12

    if(iter>1)
      nonzero_C(C+jump2*m1+m3,jump2,m1,m3);
    strassen(temp1,m2,temp2,m3,C+jump2*m1+m3,jump2,m1,m2,m3,iter,flag);   //C22+=temp1*temp2

    atomic_subtract(C+jump2*m1,jump2,C+jump2*m1+m3,jump2,m1,m3);  //C22-=C21

    //M7
    subtract(A+m2,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);   //temp1=A12-A22
    add(B+jump1*m2,jump1,B+jump1*m2+m3,jump1,temp2,m3,m2,m3);   //temp2=B21+B22

    if(iter>1)
      nonzero_C(C,jump2,m1,m3);
    strassen(temp1,m2,temp2,m3,C,jump2,m1,m2,m3,iter,0);  //C11+=temp1*temp2

    delete temp1;
    temp1=new float[m1*m3]();

    //M3
    subtract(B+m3,jump1,B+m3+jump1*m2,jump1,temp2,m3,m2,m3);  //temp2=B12-B22
    strassen(A,jump,temp2,m3,temp1,m3,m1,m2,m3,iter,1);  //temp1=A11*temp2
    atomic_add(temp1,m3,C+m3,jump2,m1,m3);  //C12=C12+temp1
    atomic_add(temp1,m3,C+jump2*m1+m3,jump2,m1,m3);  //C22=C22+temp1

    memset(temp1,0,sizeof(float)*m1*m3);
    //M4
    subtract(B+jump1*m2,jump1,B,jump1,temp2,m3,m2,m3);  //temp2=B21-B11
    strassen(A+jump*m1+m2,jump,temp2,m3,temp1,m3,m1,m2,m3,iter,1);  //temp1=A22*temp2
    atomic_add(temp1,m3,C+jump2*m1,jump2,m1,m3);  //C21=C21+temp1
    atomic_add(temp1,m3,C,jump2,m1,m3);  //C11=C11+temp1
    delete temp1,temp2;
  }
}

int nearest_ideal(int &n,int &temp)
{
  temp=0;
  int pow=1;
  int m=n;
  while(m>threshold){
    if(m%2==1){
      temp+=pow;
      m++;
    }
    m/=2;
    pow*=2;
  }
  n+=temp;
  return pow;
}

int main(int argc,char** argv)
{
  ifstream in(argv[1]);
  ifstream in1(argv[2]);
  ofstream out(argv[3]);
  float *A,*B,*C;
  int n1,n2,n3;
  in>>n1>>n2;
  in1>>n2>>n3;
  out<<n1<<"\t"<<n3<<"\n";
  int temp1,temp2,temp3;
  int iter=nearest_ideal(n1,temp1);
  iter=std::min(iter,nearest_ideal(n2,temp2));
  iter=std::min(iter,nearest_ideal(n3,temp3));

  A=new float[n1*n2];
  B=new float[n2*n3];
  C=new float[n1*n3]();

  for(int i=0; i<n1-temp1; i++)
  {
        for(int j=0; j<n2-temp2; j++)
                in>>A[i*n2+j];
	for(int j=n2-temp2;j<n2;j++)
		A[i*n2+j]=0;
  }
  for(int i=(n1-temp1)*n2;i<n1*n2;i++)
	A[i]=0;
  in.close();
  for(int i=0; i<n2-temp2; i++)
  {
        for(int j=0; j<n3-temp3; j++)
                in1>>B[i*n3+j];
	for(int j=n3-temp3;j<n3;j++)
		B[i*n3+j]=0;
  }
  for(int i=(n2-temp2)*n3;i<n2*n3;i++)
	B[i]=0;
  in1.close();

  struct timespec start,end;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start);
  strassen(A,n2,B,n3,C,n3,n1,n2,n3,iter,1);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end);

  cout<<"Time taken: "<<(end.tv_nsec-start.tv_nsec)+(1e+9)*(end.tv_sec-start.tv_sec);
  for(int i=0;i<n1-temp1;i++)
  {
    for(int j=0;j<n3-temp3;j++)
      out<<C[i*n3+j]<<"\t";
    out<<"\n";
  }
  delete C,A,B;
}
