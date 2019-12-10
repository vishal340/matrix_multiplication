#include<iostream>
#include<fstream>
#include<cmath>
#include"mpi.h"

//Please provide number of processor 
//some SMALL multiple of sqaure for better result
//like P=k*q^2 with k being small
//This code performs C+=A*B
//In this code C is defined zero before multiplication
//But you can change it to nonzero C aswell

template<typename T>
void add(T* A,int jump,T* B,int jump1,T* C,int jump2,int n1,int n2)
{
for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++)
    C[i*jump2+j]=A[i*jump+j]+B[i*jump1+j];
}
template<typename T>
void atomic_add(T* A,int jump,T* B,int jump1,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]+=A[i*jump+j];
}
template<typename T>
void subtract(T* A,int jump,T* B,int jump1,T* C,int jump2,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      C[i*jump2+j]=A[i*jump+j]-B[i*jump1+j];
}
template<typename T>
void atomic_subtract(T* A,int jump,T* B,int jump1,int n1,int n2)
{
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]-=A[i*jump+j];
}
template<typename T>
void multiply(T* A,int jump,T* B,int jump1,T* C,int jump2,int n1,int n2,int n3)
{
  for(int j=0;j<n1;j+=2){
    for(int i=0;i<n2;i+=2){
      for(int k=0;k<n3;k++){
        C[j*jump2+k]+=A[i+j*jump]*B[i*jump1+k];
        C[j*jump2+k]+=A[i+1+j*jump]*B[(i+1)*jump1+k];
        C[(j+1)*jump2+k]+=A[i+(j+1)*jump]*B[i*jump1+k];
        C[(j+1)*jump2+k]+=A[(i+1)+(j+1)*jump]*B[(i+1)*jump1+k];
      }
    }
  }
}
template<typename T>
void strassen(T* A,const int jump,T* B,const int jump1, \
              T* C,const int jump2,const int n1,const int n2,const int n3,int iter)
{
  if(iter==1)
    multiply(A,jump,B,jump1,C,jump2,n1,n2,n3);
  else
  {
    iter/=2;
    int m1=n1/2,m2=n2/2,m3=n3/2;
    T* temp1=new T[m1*m2];
    float* temp2=new T[m2*m3];
    float* temp3=new T[m1*m3]();

    add(A,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);   //temp1=A11+A22
    add(B,jump1,B+jump1*m2+m3,jump1,temp2,m3,m2,m3);   //temp2=B11+B22
    strassen(temp1,m2,temp2,m3,temp3,m3,m1,m2,m3,iter);   //temp3=temp1*temp2
    atomic_add(temp3,m3,C+jump2*m1+m3,jump2,m1,m3);  //C22+=temp3
    atomic_add(temp3,m3,C,jump2,m1,m3);  //C11+=temp3
    memset(temp3,0,sizeof(T)*m1*m3);

    add(A+jump*m1,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);  //temp1=A21+A22
    strassen(temp1,m2,B,jump1,temp3,m3,m1,m2,m3,iter);  //temp3=temp1*B11
    atomic_add(temp3,m3,C+jump2*m1,jump2,m1,m3);   //C21+=temp3
    atomic_subtract(temp3,m3,C+jump2*m1+m3,jump2,m1,m3);  //C22-=temp3
    memset(temp3,0,sizeof(T)*m1*m3);

    add(A,jump,A+m2,jump,temp1,m2,m1,m2);  //temp1=A11+A12
    strassen(temp1,m2,B+jump1*m2+m3,jump1,temp3,m3,m1,m2,m3,iter);  //temp3=temp1*B22
    atomic_add(temp3,m3,C+m3,jump2,m1,m3);    //C12+=temp3
    atomic_subtract(temp3,m3,C,jump2,m1,m3);  //C11-=temp3
    memset(temp3,0,sizeof(T)*m1*m3);

    subtract(A +jump*m1,jump,A ,jump,temp1,m2,m1,m2);   //temp1=A21-A11
    add(B,jump1,B+m3,jump1,temp2,m3,m2,m3);    //temp2=B11+B12
    strassen(temp1,m2,temp2,m3,temp3,m3,m1,m2,m3,iter);   //temp3=temp1*temp2
    atomic_add(temp3,m3,C+jump2*m1+m3,jump2,m1,m3);  //C22+=temp3
    memset(temp3,0,sizeof(T)*m1*m3);

    subtract(A+m2,jump,A+jump*m1+m2,jump,temp1,m2,m1,m2);   //temp1=A12-A22
    add(B+jump1*m2,jump1,B+jump1*m2+m3,jump1,temp2,m3,m2,m3);   //temp2=B21+B22
    strassen(temp1,m2,temp2,m3,temp3,m3,m1,m2,m3,iter);  //temp3=temp1*temp2
    atomic_add(temp3,m3,C,jump2,m1,m3);   //C11+=temp3
    memset(temp3,0,sizeof(T)*m1*m3);
    delete temp1;

    subtract(B+m3,jump1,B+m3+jump1*m2,jump1,temp2,m3,m2,m3);  //temp2=B12-B22
    strassen(A,jump,temp2,m3,temp3,m3,m1,m2,m3,iter);  //temp3=A11*temp2
    atomic_add(temp3,m3,C+m3,jump2,m1,m3);  //C12+=temp3
    atomic_add(temp3,m3,C+jump2*m1+m3,jump2,m1,m3);  //C22+=temp3
    memset(temp3,0,sizeof(T)*m1*m3);

    subtract(B+jump1*m2,jump1,B,jump1,temp2,m3,m2,m3);  //temp2=B21-B11
    strassen(A+jump*m1+m2,jump,temp2,m3,temp3,m3,m1,m2,m3,iter);  //temp3=A22*temp2
    atomic_add(temp3,m3,C+jump2*m1,jump2,m1,m3);  //C21+=temp3
    atomic_add(temp3,m3,C,jump2,m1,m3);  //C11+=temp3
    delete temp2,temp3;
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
  n+=temp;
  temp+=t;
  if(m%2==1){
    n+=(pow*temp1);
    temp+=(pow*temp1);
  }
  return pow;
}

int main(int argc,char** argv){
  const int threshold=32;   //when the iterative method of strassen reaches this size
                              //It performs normal multiplication on the small matrix
  MPI_Init(&argc,&argv);
  int rank,NUM_processors;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&NUM_processors);
  int iter[3]={NUM_processors,1,1};
  int temp[3]={0,0,0};
  int n[3]; //the three dimensions whole file and per processor respectively
  int row,cont,col;  //This is useful dimension when we have to append zeros
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
    float quotient=std::sqrt((n[1]*n[1])/(n[0]*n[2]));
    float difference=NUM_processors*quotient-1;
    iter[0]=NUM_processors,iter[1]=1;
    for(int i=2;i<=(int)(std::sqrt(NUM_processors));i++){
      if(NUM_processors%i == 0){
        if(std::fabs((NUM_processors/i)*quotient-i)<difference){
          iter[0]=(int)(NUM_processors/i);
          iter[1]=i;
        }
      }
    }
    iter[2]=nearest_ideal(n[0],temp[0],iter[0],threshold);
    iter[2]=std::min(iter[2],nearest_ideal(n[1],temp[1],iter[1],threshold));
    iter[2]=std::min(iter[2],nearest_ideal(n[2],temp[2],iter[0],threshold));
  }
  MPI_Bcast(n,3,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(iter,3,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(temp,3,MPI_INT,0,MPI_COMM_WORLD);
  if(rank==(NUM_processors-1)){
    std::fstream out;
    out.open(argv[3],std::ios::out|std::ios::binary);
    out.write((char*)&n[0],sizeof(int));
    out.write((char*)&n[2],sizeof(int));
    int fill=0;
    out.seekp(sizeof(float)*(n[0]*n[2]-1),std::ios::cur);
    out.write((char*)&fill,sizeof(int));
    out.close();
  }
  const int dim_rank[2]={(int)(rank/iter[1]), rank%iter[1]};
  const int m[3]={(int)(n[0]/iter[0]),(int)(n[1]/iter[1]),(int)(n[2]/iter[0])};
  n[0]-=temp[0];n[1]-=temp[1];n[2]-=temp[2];
  const int gap[3]={n[0]%iter[0],n[1]%iter[1],n[2]%iter[0]};
  const int m1[3]={(int)(n[0]/iter[0]),(int)(n[1]/iter[1]),(int)(n[2]/iter[0])};

  MPI_Status status;
  MPI_Request request;

  float* matrix_loc1=new float[m[0]*m[1]]();  //matrix A for local processor
  float* matrix_loc2=new float[m[1]*m[2]]();  //matrix B for local processor
    
  if(temp[0] !=0 || temp[1] !=0 || temp[2] !=0){
    row=m1[0]+(gap[0] > dim_rank[0] ? 1:0); //no. of rows for output matrix from this processor
    cont=m1[1]+((gap[1]>dim_rank[1]) ? 1:0);
    int offset1=((gap[0] > dim_rank[0] ? dim_rank[0]:gap[0])+dim_rank[0]*m1[0])*n[1]+\
                ((gap[1]>dim_rank[1]) ? dim_rank[1]:gap[1])+dim_rank[1]*m1[1];
    MPI_File matrix1;
    MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix1);
    for(int i=0;i<row;i++){
      MPI_File_iread_at(matrix1,2*sizeof(int)+sizeof(float)*(i*n[1]+offset1),matrix_loc1+i*m[1],cont,MPI_FLOAT,&request);
      MPI_Wait(&request,&status);
    }
    MPI_File_close(&matrix1);
    
    col=m1[2]+(gap[2] > (dim_rank[0]+dim_rank[1])%iter[0] ? 1:0); //no. of col for output matrix from this processor
    int offset2=(gap[1]>dim_rank[1] ? (dim_rank[1]*(m1[1]+1)*n[2]):((gap[1]+dim_rank[1]*m1[1])*n[2]))+\
                (gap[2]>(dim_rank[0]+dim_rank[1])%iter[0] ? (((dim_rank[0]+dim_rank[1])%iter[0])*(m1[2]+1)):\
                (gap[2]+((dim_rank[0]+dim_rank[1])%dim_rank[0])*m1[2]));
    MPI_File matrix2;
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix2);
    for(int i=0;i<cont;i++){
      MPI_File_iread_at(matrix2,2*sizeof(int)+sizeof(float)*(i*n[2]+offset2),matrix_loc2+i*m[2],col,MPI_FLOAT,&request);
      MPI_Wait(&request,&status);
    }
    MPI_File_close(&matrix2);
  }
  else{
    MPI_File matrix1;
    MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix1);
    MPI_Datatype new_float1;
    MPI_Type_vector(m[0],m[1],n[1],MPI_FLOAT,&new_float1);
    MPI_Type_commit(&new_float1);
    MPI_File_set_view(matrix1,2*sizeof(int)+sizeof(float)*(dim_rank[0]*n[1]*m[0]+ \
                    dim_rank[1]*m[1]),MPI_FLOAT,new_float1,"native",MPI_INFO_NULL);
    MPI_File_iread(matrix1,matrix_loc1,m[0]*m[1],MPI_FLOAT,&request);
    MPI_Wait(&request,&status);
    MPI_File_close(&matrix1);

    MPI_File matrix2;
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix2);
    MPI_Datatype new_float2;
    MPI_Type_vector(m[1],m[2],n[2],MPI_FLOAT,&new_float2);
    MPI_Type_commit(&new_float2);
    MPI_File_set_view(matrix2,2*sizeof(int)+sizeof(float)*(dim_rank[1]*n[2]*m[1]+((dim_rank[0]+dim_rank[1])%iter[0])*m[2]),\
                    MPI_FLOAT,new_float2,"native",MPI_INFO_NULL);
    MPI_File_iread(matrix2,matrix_loc2,m[1]*m[2],MPI_FLOAT,&request);
    MPI_Wait(&request,&status);
    MPI_File_close(&matrix2);
  }

  float** matrix_loc3;      //matrix C for local processor
  double start,end;

  MPI_Barrier(MPI_COMM_WORLD);
  if(iter[1]>1){
    const int iteration[2]={(int)(iter[0]/iter[1]),iter[0]%iter[1]};
    matrix_loc3=new float*[iteration[0]+(iteration[1]>0 ? 1:0)];
    for(int i=0;i<iteration[0]+(iteration[1]>0 ? 1:0);i++)
      matrix_loc3[i]=new float[m[0]*m[2]]();
    start=MPI_Wtime();
    const int source1=((dim_rank[0]+iter[0]-1)%iter[0])*iter[1]+(dim_rank[1]+1)%iter[1];
    const int dest1=((dim_rank[0]+iter[0]-iter[1]+1)%iter[0])*iter[1]+(dim_rank[1]+iter[1]-1)%iter[1];
    const int source2=(dim_rank[1]+1)%iter[1]+dim_rank[0]*iter[1];
    const int dest2=(dim_rank[1]+iter[1]-1)%iter[1]+dim_rank[0]*iter[1];
    const int source3=((dim_rank[0]+iter[1])%iter[0])*iter[1]+dim_rank[1];
    const int dest3=((dim_rank[0]+iter[0]-iter[1])%iter[0])*iter[1]+dim_rank[1];
    for(int i=0;i<(iteration[0]-3);i+=2){
      for(int j=0;j<(iter[1]-1);j++){
        strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i],m[2],m[0],m[1],m[2],iter[2]);
        MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest1,0,source1,0,MPI_COMM_WORLD,&status);
        MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,dest2,1,source2,1,MPI_COMM_WORLD,&status);
      }
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i],m[2],m[0],m[1],m[2],iter[2]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      for(int j=0;j<(iter[1]-1);j++){
        strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i+1],m[2],m[0],m[1],m[2],iter[2]);
        MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,source1,0,dest1,0,MPI_COMM_WORLD,&status);
        MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,source2,1,dest2,1,MPI_COMM_WORLD,&status);
      }
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i+1],m[2],m[0],m[1],m[2],iter[2]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
    }
    if(iteration[0]>1){
      int i;
      if(iteration[0]%2==0){i=iteration[0]-2;}
      else{i=iteration[0]-3;}
    for(int j=0;j<(iter[1]-1);j++){
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i],m[2],m[0],m[1],m[2],iter[2]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest1,0,source1,0,MPI_COMM_WORLD,&status);
      MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,dest2,1,source2,1,MPI_COMM_WORLD,&status);
    }
    strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i],m[2],m[0],m[1],m[2],iter[2]);
    MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
    for(int j=0;j<(iter[1]-1);j++){
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i+1],m[2],m[0],m[1],m[2],iter[2]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,source1,0,dest1,0,MPI_COMM_WORLD,&status);
      MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,source2,1,dest2,1,MPI_COMM_WORLD,&status);
    }
    strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i+1],m[2],m[0],m[1],m[2],iter[2]);
    }
    if(iteration[0]%2==1){
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      for(int j=0;j<(iter[1]-1);j++){
        strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[iteration[0]-1],m[2],m[0],m[1],m[2],iter[2]);
        MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest1,0,source1,0,MPI_COMM_WORLD,&status);
        MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,dest2,1,source2,1,MPI_COMM_WORLD,&status);
      }
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[iteration[0]-1],m[2],m[0],m[1],m[2],iter[2]);
    }
    if(iteration[1]>0){
      const int source6=(dim_rank[1]<iteration[1]? source3:((dim_rank[0]+iter[1]-\
                        iteration[1])%iter[0])*iter[1]+dim_rank[1]);
      const int dest6=(dim_rank[1]<iteration[1]? dest3:((dim_rank[0]+iter[0]\
                      -iter[1]+iteration[1])%iter[0])*iter[1]+dim_rank[1]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest6,2,source6,2,MPI_COMM_WORLD,&status);

      const int source4=dim_rank[1]<iteration[1]? ((dim_rank[1]+1)%iteration[1]+\
                        dim_rank[0]*iter[1]):(dim_rank[0]*iter[1]+dim_rank[1]);
      const int dest4=dim_rank[1]<iteration[1]? (dim_rank[1]+iteration[1]-1)%iteration[1]+\
                      dim_rank[0]*iter[1]:(dim_rank[0]*iter[1]+dim_rank[1]);
      const int source5=dim_rank[1]<iteration[1]? ((dim_rank[1]+iteration[1]-1)%iteration[1]+((dim_rank[0]+\
                        iter[0]-1)%iter[0])*iter[1]):(((dim_rank[0]+1)%iter[0])*iter[1]+dim_rank[1]);
      const int dest5=dim_rank[1]<iteration[1]? (dim_rank[1]+1)%iteration[1]+((dim_rank[0]+1)%iter[0])\
                      *iter[1]:(((dim_rank[0]+iter[0]-1)%iter[0])*iter[1]+dim_rank[1]);
      for(int i=0;i<iteration[1];i++){
          strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[iteration[0]],m[2],m[0],m[1],m[2],iter[2]);
          MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest5,0,source5,0,MPI_COMM_WORLD,&status);

          MPI_Sendrecv_replace(matrix_loc1,m[0]*m[1],MPI_FLOAT,dest4,1,source4,1,MPI_COMM_WORLD,&status);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end=MPI_Wtime();
    if(temp[0]==0 && temp[1]==0 && temp[2]==0){
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      MPI_Datatype new_float3;
      MPI_Type_vector(m[0],m[2],n[2],MPI_FLOAT,&new_float3);
      MPI_Type_commit(&new_float3);
      for(int i=0;i<iteration[0];i++){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(((dim_rank[0]+dim_rank[1]+i*iter[1])%iter[0])\
                          *m[2]+dim_rank[1]*n[2]*m[0]),MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_iwrite(matrix3,matrix_loc3[i],m[0]*m[2],MPI_FLOAT,&request);
        MPI_Wait(&request,&status);
      }
      if(dim_rank[1]<iteration[1]){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(((dim_rank[0]+dim_rank[1]+iteration[0]*iter[1])%iter[0])\
                          *m[2]+dim_rank[1]*n[2]*m[0]),MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_iwrite(matrix3,matrix_loc3[iteration[0]],m[0]*m[2],MPI_FLOAT,&request);
        MPI_Wait(&request,&status);
      }
      MPI_File_close(&matrix3);
    }
    else{
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      for(int i=0;i<iteration[0];i++){
        int shift=((dim_rank[0]+dim_rank[1]+i*iter[1])%iter[0]);
        for(int j=0;j<row;j++){
          MPI_File_iwrite_at(matrix3,2*sizeof(int)+sizeof(float)*((gap[2]>shift?shift*(m1[2]+1):(gap[2]+shift*m1[2]))+\
                            (gap[1]>dim_rank[1]? dim_rank[1]*(m1[1]+1)*n[2]:(gap[1]+dim_rank[1])*n[2])),
                            matrix_loc3[i]+j*m[2],col,MPI_FLOAT,&request);
          MPI_Wait(&request,&status);
      }
    }
    if(dim_rank[1]<iteration[1]){
      int shift=((dim_rank[0]+dim_rank[1]+iteration[0]*iter[1])%iter[0]);
      for(int j=0;j<row;j++){
        MPI_File_iwrite_at(matrix3,2*sizeof(int)+sizeof(float)*((gap[2]>shift?shift*(m1[2]+1):(gap[2]+shift*m1[2]))+\
                          (gap[1]>dim_rank[1]? dim_rank[1]*(m1[1]+1)*n[2]:(gap[1]+dim_rank[1])*n[2])),
                          matrix_loc3[iteration[0]]+j*m[2],col,MPI_FLOAT,&request);
        MPI_Wait(&request,&status);
      }
    }
    MPI_File_close(&matrix3);
    }
  }
  else
  {
    matrix_loc3=new float*[NUM_processors];
    for(int i=0;i<NUM_processors;i++)
      matrix_loc3[i]=new float[m[0]*m[2]]();
    start=MPI_Wtime();
    const int source3=(dim_rank[0]+1)%NUM_processors;
    const int dest3=(dim_rank[0]+NUM_processors-1)%NUM_processors;
    for(int i=0;i<(NUM_processors-1);i++){
      strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[i],m[2],m[0],m[1],m[2],iter[2]);
      MPI_Sendrecv_replace(matrix_loc2,m[1]*m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
    }
    strassen(matrix_loc1,m[1],matrix_loc2,m[2],matrix_loc3[NUM_processors-1],m[2],m[0],m[1],m[2],iter[2]);
    MPI_Barrier(MPI_COMM_WORLD);
    end=MPI_Wtime();
    if(temp[0]==0 && temp[1]==0 && temp[2]==0){
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      MPI_Datatype new_float3;
      MPI_Type_vector(m[0],m[2],n[2],MPI_FLOAT,&new_float3);
      MPI_Type_commit(&new_float3);
      for(int i=0;i<NUM_processors;i++){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(dim_rank[0]*n[2]*m[0]+((i+dim_rank[0])%NUM_processors)*m[2]),\
                          MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_iwrite(matrix3,matrix_loc3[i],m[0]*m[2],MPI_FLOAT,&request);
        MPI_Wait(&request,&status);
      }
      MPI_File_close(&matrix3);
    }
    else{
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      int offset3=((gap[0] > dim_rank[0] ? dim_rank[0]:gap[0])+dim_rank[0]*m1[0])*n[2];
      for(int i=0;i<NUM_processors;i++){
        int k=(i+dim_rank[0])%NUM_processors;
        int offset_loop=offset3+(gap[2]>k ? k*(m1[2]+1):gap[2]+k*m1[2]);
        for(int j=0;j<row;j++){
          MPI_File_iwrite_at(matrix3,2*sizeof(int)+sizeof(float)*(offset_loop+j*n[2]),\
                            matrix_loc3[i]+j*m[2],col,MPI_FLOAT,&request);
          MPI_Wait(&request,&status);
      }
    }
    MPI_File_close(&matrix3);
    }
  }
  
  std::cout<<matrix_loc3[0][0]<<'\t'<<matrix_loc3[0][1]<<'\t'<<matrix_loc3[0][2]<<'\t';
  if(rank==0){
    std::cout<<"Time taken: "<<end-start<<'\n';
  }
  delete matrix_loc1,matrix_loc2,matrix_loc3;
  MPI_Finalize();
  return 0;
}
