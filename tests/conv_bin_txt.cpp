#include<iostream>
#include<fstream>
using namespace std;
int main(int argc,char** argv)
{
    ifstream in(argv[1],ios::binary);
    ofstream out(argv[2]);
    int n,m;
    in.read((char*)&n,sizeof(int));
    in.read((char*)&m,sizeof(int));
    out<<n<<"\t"<<m<<"\n";
    float *a=new float[n*m];
    in.read((char*)a,sizeof(float)*n*m);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<m;j++)
            out<<a[i*m+j]<<"\t";
        out<<"\n";
    }
    delete a;
}
