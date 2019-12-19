#include<iostream>
#include<fstream>

int main(int argc,char** argv){
    std::fstream in,out;
    in.open(argv[1],std::ios::in|std::ios::binary);
    out.open(argv[2],std::ios::out);
    int n1,n2;float f;
    in.read((char*)&n1,sizeof(int));
    out<<n1<<'\t';
    in.read((char*)&n2,sizeof(int));
    out<<n2<<'\n';
    for(int j=0;j<n1;j++){
        for(int i=0;i<n2;i++){
            in.read((char*)&f,sizeof(float));
            out<<f<<'\t';
        }
        out<<'\n';
    }
    return 0;
}
