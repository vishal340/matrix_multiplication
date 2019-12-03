#include<iostream>
#include<fstream>
#include<cmath>
#include"mpi.h"

//Please provide number of processor 
//some SMALL multiple of sqaure for better result
//like P=k*q^2 with k being small

void add(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2)
{
for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++)
    C[i*jump2+j]=A[i*jump+j]+B[i*jump1+j];
}

void atomic_add(float* A,int jump,float* B,int jump1,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]+=A[i*jump+j];
}

void subtract(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      C[i*jump2+j]=A[i*jump+j]-B[i*jump1+j];
}

void atomic_subtract(float* A,int jump,float* B,int jump1,int n1,int n2)
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

void nonzero_C(float* C,int jump2,int m1,int m3)
{
  atomic_add(C+m3/2,jump2,C,jump2,m1/2,m3/2);
  atomic_add(C+jump2*m1/2,jump2,C+jump2*m1/2+m3/2,jump2,m1/2,m3/2);
  atomic_subtract(C,jump2,C+jump2*m1/2+m3/2,jump2,m1/2,m3/2);
}

void strassen(float* A,const int jump,float* B,const int jump1, \
              float* C,const int jump2,const int n1,const int n2,const int n3,int iter,const bool flag)
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

int nearest_ideal(int &n,int &temp,const int temp1,const int threshold)
{
    int t=(temp1-n%temp1)%temp1;
    int pow=1;
    n+=t;
    int m=(int)(n/temp1);
    while(m>threshold){
        if(m%2==1){
            temp+=pow;
            m++;
        }
    m/=2;
    pow*=2;
    }
    temp*=temp1;
    temp+=t;
    n+=temp;
    return pow;
}

int main(int argc,char** argv){
    const int threshold=32;   //when the iterative method of strassen reaaches this size
                              //It performs normal multiplication on the small matrix
    MPI_Init(&argc,&argv);
    int rank,NUM_processors;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&NUM_processors);
    int* iter=new int[3];
    int* temp=new int[3]();
    int* n=new int[3]; //the three dimensions whole file and per processor respectively
    if(NUM_processors<2)
        throw("number of processors should be more than 1");
    if(rank==0){
        std::fstream in;
        in.open(argv[1],std::ios::in|std::ios::binary);
        in.read((char*)&n[0],sizeof(int));
        in.seekg(sizeof(int),std::ios::beg);
        in.read((char*)&n[1],sizeof(int));
        in.close();
        std::fstream in1;
        in1.open(argv[2],std::ios::in|std::ios::binary);
        in1.read((char*)&n[1],sizeof(int));
        in1.seekg(sizeof(int),std::ios::beg);
        in1.read((char*)&n[2],sizeof(int));
        in1.close();
        iter[0]=NUM_processors,iter[1]=1;
        for(int i=2;i <= (int)(std::sqrt(NUM_processors));i++){
              if(NUM_processors%(i*i) == 0){
                  iter[0]=NUM_processors/(i*i);
                  iter[1]=i;
              }
        }
        iter[2]=nearest_ideal(n[0],temp[0],iter[0]*iter[1],threshold);
        iter[2]=std::min(iter[2],nearest_ideal(n[1],temp[1],iter[1],threshold));
        iter[2]=std::min(iter[2],nearest_ideal(n[2],temp[2],iter[0]*iter[1],threshold));
    }
    //std::cout<<n[0]<<'\t'<<n[1]<<'\t'<<n[2]<<'\n';
    MPI_Bcast(n,3,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(iter,3,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(temp,3,MPI_INT,0,MPI_COMM_WORLD);
    //std::cout<<iter[0]<<'\t'<<iter[1]<<'\t'<<iter[2]<<'\n';
    const int dim_rank[2]={(int)(rank/iter[1]), rank%iter[1]};
    const int m[3]={(int)(n[0]/(iter[0]*iter[1])),(int)(n[1]/iter[1]),(int)(n[2]/(iter[0]*iter[1]))};

    float* matrix_loc1=new float[m[0]*m[1]];
    MPI_File matrix1;
    MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix1);
    MPI_Datatype new_float1;
    MPI_Type_vector(m[0],m[1],n[1],MPI_FLOAT,&new_float1);
    MPI_Type_commit(&new_float1);
    MPI_File_set_view(matrix1,2*sizeof(int)+sizeof(float)*(dim_rank[0]*n[1]*m[0]+dim_rank[1]*m[1]),MPI_FLOAT,new_float1,"native",MPI_INFO_NULL);
    MPI_File_read_all(matrix1,matrix_loc1,m[0]*m[1],MPI_FLOAT,MPI_STATUS_IGNORE);
    MPI_File_close(&matrix1);

    float* matrix_loc2=new float[m[1]*m[2]];
    MPI_File matrix2;
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix2);
    MPI_Datatype new_float2;
    MPI_Type_vector(m[1],m[2],n[2],MPI_FLOAT,&new_float2);
    MPI_Type_commit(&new_float2);
    MPI_File_set_view(matrix2,2*sizeof(int)+sizeof(float)*((dim_rank[0]%iter[0])*m[2]*iter[1]+dim_rank[1]*m[1]*n[2] \
                      +((dim_rank[0]%iter[0]+dim_rank[1])%iter[1])*m[2]),MPI_FLOAT,new_float2,"native",MPI_INFO_NULL);
    MPI_File_read_all(matrix2,matrix_loc2,m[1]*m[2],MPI_FLOAT,MPI_STATUS_IGNORE);
    MPI_File_close(&matrix2);

    float** matrix_loc3=new float*[iter[0]];
    for(int i=0;i<iter[0];i++)
    matrix_loc3[i]=new float[m[0]*m[2]]();

    double start=MPI_Wtime();
    for(int i=0;i<iter[0];i++)
    {
      int k=(i+(int)(dim_rank[0]/iter[1]))%iter[0];
      for(int j=0;j<iter[1];j++){
        MPI_Status status1,status2;
          strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[k],m[2],m[0],m[1],m[2],iter[2],1);
          MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,(rank+iter[1]-1)%iter[1]+dim_rank[0]*iter[1],0, \
                        (rank+1)%iter[1]+dim_rank[0]*iter[1],0,MPI_COMM_WORLD,&status1);
          MPI_Sendrecv_replace(matrix_loc3[k],m[0]*m[2],MPI_FLOAT,(rank+iter[1]-1)%iter[1]+dim_rank[0]*iter[1],0, \
                        (rank+1)%iter[1]+dim_rank[0]*iter[1],0,MPI_COMM_WORLD,&status2);
      }
      MPI_Status status3;
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,(rank+NUM_processors-iter[1]*iter[1])%NUM_processors,0, \
                    (rank+iter[1]*iter[1])%NUM_processors,0,MPI_COMM_WORLD,&status3);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double end=MPI_Wtime();

    if(rank==0){
        std::cout<<"Time taken: "<<end-start<<'\n';
        std::fstream out;
        out.open(argv[3],std::ios::out|std::ios::binary);
        out.write((char*)&n[0],sizeof(int));
        out.seekp(sizeof(int),std::ios::beg);
        out.write((char*)&n[2],sizeof(int));
        int fill=0;
        out.seekp(sizeof(float)*(m[0]*m[2]-1),std::ios::cur);
        out.write((char*)&fill,sizeof(int));
        out.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File matrix3;
    MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
    MPI_Datatype new_float3;
    MPI_Type_vector(m[0],m[2],n[2],MPI_FLOAT,&new_float3);
    MPI_Type_commit(&new_float3);
    for(int i=0;i<iter[0];i++){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(m[2]*iter[1]*i+dim_rank[0]*n[2]*m[0]+ \
                          ((dim_rank[0]%iter[0]+dim_rank[1])%iter[1])*m[2]),MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_write_all(matrix3,matrix_loc3[i],m[0]*m[2],MPI_FLOAT,MPI_STATUS_IGNORE);
    }
    MPI_File_close(&matrix3);

    MPI_Finalize();
    return 0;
}
