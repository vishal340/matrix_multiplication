#include<iostream>
#include<fstream>

using namespace std;

int main(int argc,char** argv)
{
    ifstream in(argv[1]);
    ifstream in1(argv[2]);
    ofstream out(argv[3]);
    int n1,n2,n3;
    in>>n1>>n2;
    in1>>n2>>n3;
    out<<n1<<"\t"<<n3<<"\n";
    float *A=new float[n1*n2];
    float *B=new float[n2*n3];
    float *C=new float[n1*n3]();
    for(int i=0;i<n1;i++)
        for(int j=0;j<n2;j++)
            in>>A[i*n2+j];
    for(int i=0;i<n2;i++)
        for(int j=0;j<n3;j++)
            in1>>B[i*n3+j];
    in.close();
    in1.close();
    for(int i=0;i<n1;i++)
        for(int k=0;k<n2;k++)
            for(int j=0;j<n3;j++)
                C[i*n3+j]+=A[i*n2+k]*B[k*n3+j];
    for(int i=0;i<n1;i++)
    {
        for(int j=0;j<n3;j++)
            out<<C[i*n3+j]<<"\t";
        out<<"\n";
    }
    delete A,B,C;
}
