#include<iostream>
#include<fstream>
//#include"mpi.h"

int main(int argc,char** argv){
    std::fstream in1,in2;
    float f1,f2;
    int n1,n2;
    in1.open(argv[1],std::ios::in);
    in2.open(argv[2],std::ios::in);
    in1>>n1;
    in1>>n2;
    int count=0;
    for(int i=0;i<n1;i++){
        for(int j=0;j<n2;j++){
            in1>>f1;
            in2>>f2;
            if(std::abs(f1-f2)>10){
                count++;
                //std::cout<<i<<" "<<j<<'\t';
            }
        }
    }
    std::cout<<count<<'\n';
    return 0;
}
