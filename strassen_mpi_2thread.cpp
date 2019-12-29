#include<iostream>
#include<fstream>
#include<cmath>
#include<memory>
#include<thread>
#include<mpi.h>
#include<cstring>

using std::unique_ptr;
using std::thread;

//Please provide number of processor 
//some SMALL multiple of sqaure for better result
//like P=k*q^2 with k being small
//This code performs C+=A*B
//In this code C is defined zero before multiplication
//But you can change it to nonzero C aswell

template<typename T>
void add(const T* const A,int jump,const T* const B,int jump1,T* C,int jump2,int n1,int n2){
for(int i=0;i<n1;i++)
  for(int j=0;j<n2;j++)
    C[i*jump2+j]=A[i*jump+j]+B[i*jump1+j];
}
template<typename T>
void atomic_add(const T* const A,int jump,T* B,int jump1,int n1,int n2){
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]+=A[i*jump+j];
}
template<typename T>
void subtract(const T* const A,int jump,const T* const B,int jump1,T* C,int jump2,int n1,int n2){
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      C[i*jump2+j]=A[i*jump+j]-B[i*jump1+j];
}
template<typename T>
void atomic_subtract(const T* const A,int jump,T* B,int jump1,int n1,int n2){
  for(int i=0;i<n1;i++)
    for(int j=0;j<n2;j++)
      B[i*jump1+j]-=A[i*jump+j];
}
template<typename T>
void multiply(const T* const A,int jump,const T* const B,int jump1,T* C,int jump2,int n1,int n2,int n3){
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
void C_adjust1(T* C,int jump,int n1,int n2){
  for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++){
      C[i*jump+j]+=(C[i*jump+j+n2]-C[(i+n1)*jump+j]);
      C[(i+n1)*jump+j+n2]-=C[i*jump+j];
    }
  }
}
template<typename T>
void C_adjust2(T* C,int jump,int n1,int n2){
  for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++){
      C[(i+n1)*jump+j+n2]+=(C[i*jump+j]-C[(i+n1)*jump+j]+C[i*jump+j+n2]);
    }
  }
}

template<typename T>
void strassen(const T* const A,const int jump,const T* const B,const int jump1, \
              T* C,const int jump2,int m1,int m2,int m3,int iter){
  if(iter==1)
    multiply(A,jump,B,jump1,C,jump2,m1,m2,m3);
  else{
    iter>>=1;
    m1>>=1,m2>>=1,m3>>=1;
    unique_ptr<T[]>temp1=std::make_unique<T[]>(m1*m2);
    unique_ptr<T[]>temp2=std::make_unique<T[]>(m2*m3);

    C_adjust1(C,jump2,m1,m3);    //C11+=C12-C21
                                //C22-=C11
     //M5
    add(&A[0],jump,&A[m2],jump,&temp1[0],m2,m1,m2);  //temp1=A11+A12
    strassen(&temp1[0],m2,&B[jump1*m2+m3],jump1,&C[m3],jump2,m1,m2,m3,iter);  //C12+=temp1*B22
    atomic_subtract(&C[m3],jump2,&C[0],jump2,m1,m3);  //C11-=C12
    //M3
    subtract(&B[m3],jump1,&B[m3+jump1*m2],jump1,&temp2[0],m3,m2,m3);  //temp2=B12-B22
    strassen(&A[0],jump,&temp2[0],m3,&C[m3],jump2,m1,m2,m3,iter);  //C12+=A11*temp2
    //M4
    subtract(&B[jump1*m2],jump1,&B[0],jump1,&temp2[0],m3,m2,m3);  //temp2=B21-B11
    strassen(&A[jump*m1+m2],jump,&temp2[0],m3,&C[jump2*m1],jump2,m1,m2,m3,iter);  //C21+=A22*temp2
    atomic_add(&C[jump2*m1],jump2,&C[0],jump2,m1,m3);   //C11+=C21
    //M2
    add(&A[jump*m1],jump,&A[jump*m1+m2],jump,&temp1[0],m2,m1,m2);  //temp1=A21+A22
    strassen(&temp1[0],m2,&B[0],jump1,&C[jump2*m1],jump2,m1,m2,m3,iter);  //C21+=temp1*B11
    //M1
    add(&A[0],jump,&A[jump*m1+m2],jump,&temp1[0],m2,m1,m2);   //temp1=A11+A22
    add(&B[0],jump1,&B[jump1*m2+m3],jump1,&temp2[0],m3,m2,m3);   //temp2=B11+B22
    strassen(&temp1[0],m2,&temp2[0],m3,&C[0],jump2,m1,m2,m3,iter);   //C11+=temp1*temp2

    C_adjust2(&C[0],jump2,m1,m3);     //C22+=C11-C21+C12
    
    //M6
    subtract(&A[jump*m1],jump,&A[0] ,jump,&temp1[0],m2,m1,m2);   //temp1=A21-A11
    add(&B[0],jump1,&B[m3],jump1,&temp2[0],m3,m2,m3);    //temp2=B11+B12
    strassen(&temp1[0],m2,&temp2[0],m3,&C[jump2*m1+m3],jump2,m1,m2,m3,iter);  //C22+=M6
    //M7
    subtract(&A[m2],jump,&A[jump*m1+m2],jump,&temp1[0],m2,m1,m2);   //temp1=A12-A22
    add(&B[jump1*m2],jump1,&B[jump1*m2+m3],jump1,&temp2[0],m3,m2,m3);   //temp2=B21+B22
    strassen(&temp1[0],m2,&temp2[0],m3,&C[0],jump2,m1,m2,m3,iter);  //C11+=M7
  }
}
int nearest_ideal(int &n,int &temp,const int temp1,int threshold)
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
  m>>=1;
  pow<<=1;
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
  int provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
  int rank,NUM_processors,flag;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&NUM_processors);
  int iter[3]={NUM_processors,1,1};
  int temp[3]={0,0,0};
  int n[3];      //the three dimensions whole file and per processor respectively
  int row,cont;  //This is useful dimension when we have to append zeros
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
    float ratio=std::sqrt((n[1]*n[1])/(n[0]*n[2]));
    float difference=NUM_processors*ratio-1;
    iter[0]=NUM_processors,iter[1]=1;
    for(int i=2;i<=(int)(std::sqrt(NUM_processors));i++){
      if(NUM_processors%i == 0){
        if(std::fabs((NUM_processors/i)*ratio-i)<difference){
          iter[0]=(int)(NUM_processors/i);
          iter[1]=i;
        }
      }
    }
    int threshold=(int)std::pow(2,(int)(-2-std::log2((double)(8*iter[1])/(double)n[1]+(double)(5*iter[0])/(double)(n[0]+n[2]))));
    if(threshold<1){
      threshold=1;
    }
    //std::cout<<"threshold: "<<threshold<<'\n';
    iter[2]=nearest_ideal(n[0],temp[0],iter[0],(n[0]/(4*iter[0]*threshold)));
    if(iter[2]>1){
      iter[2]>>=1;
    }
    iter[2]=std::min(iter[2],nearest_ideal(n[2],temp[2],iter[0],(n[2]/(4*iter[0]*threshold))));
    iter[2]=std::min(iter[2],nearest_ideal(n[1],temp[1],iter[1],(n[1]/(4*iter[1]*threshold))));
    //std::cout<<"iterations: "<<iter[2]<<'\n';
  }
  MPI_Bcast(n,3,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(iter,3,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Bcast(temp,3,MPI_INT,0,MPI_COMM_WORLD);
  
  const int dim_rank[2]={(int)(rank/iter[1]), rank%iter[1]};
  const int m[3]={(int)(n[0]/iter[0]),(int)(n[1]/iter[1]),(int)(n[2]/iter[0])};
  n[0]-=temp[0];n[1]-=temp[1];n[2]-=temp[2];
  const int gap[3]={n[0]%iter[0],n[1]%iter[1],n[2]%iter[0]};
  const int m1[3]={(int)(n[0]/iter[0]),(int)(n[1]/iter[1]),(int)(n[2]/iter[0])};

  int temp_m0=m[0]>>1;

  /*if(rank==(NUM_processors-1)){
    std::fstream out;
    out.open(argv[3],std::ios::out|std::ios::binary);
    out.write((char*)&n[0],sizeof(int));
    out.write((char*)&n[2],sizeof(int));
    int fill=0;
    out.seekp(sizeof(float)*(n[0]*n[2]-1),std::ios::cur);
    out.write((char*)&fill,sizeof(int));
    out.close();
  }*/

  MPI_Status status;
  MPI_Request request;

  unique_ptr<float[]> matrix_loc1=std::make_unique<float[]>(m[0]*m[1]);  //matrix A for local processor
  unique_ptr<float[]> matrix_loc2=std::make_unique<float[]>(m[1]*m[2]);  //matrix B for local processor
    
  if(temp[0] !=0 || temp[1] !=0 || temp[2] !=0){
    row=m1[0]+(gap[0] > dim_rank[0] ? 1:0); //no. of rows for output matrix from this processor
    cont=m1[1]+((gap[1]>dim_rank[1]) ? 1:0);
    int offset1=((gap[0] > dim_rank[0] ? dim_rank[0]:gap[0])+dim_rank[0]*m1[0])*n[1]+\
                ((gap[1]>dim_rank[1]) ? dim_rank[1]:gap[1])+dim_rank[1]*m1[1];
    MPI_File matrix1;
    MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix1);
    for(int i=0;i<row;i++){
      MPI_File_iread_at(matrix1,2*sizeof(int)+sizeof(float)*(i*n[1]+offset1),&matrix_loc1[i*m[1]],cont,MPI_FLOAT,&request);
      MPI_Wait(&request,&status);
    }
    MPI_File_close(&matrix1);
    
    int offset2=(gap[1]>dim_rank[1] ? (dim_rank[1]*(m1[1]+1)*n[2]):((gap[1]+dim_rank[1]*m1[1])*n[2]))+\
                (gap[2]>(dim_rank[0]+dim_rank[1])%iter[0] ? (((dim_rank[0]+dim_rank[1])%iter[0])*(m1[2]+1)):\
                (gap[2]+((dim_rank[0]+dim_rank[1])%iter[0])*m1[2]));
    MPI_File matrix2;
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix2);
    for(int i=0;i<cont;i++){
      MPI_File_iread_at(matrix2,2*sizeof(int)+sizeof(float)*(i*n[2]+offset2),&matrix_loc2[i*m[2]],m1[2]+(gap[2] >\
                         (dim_rank[0]+dim_rank[1])%iter[0] ? 1:0),MPI_FLOAT,&request);
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
    MPI_File_iread(matrix1,&matrix_loc1[0],m[0]*m[1],MPI_FLOAT,&request);
    MPI_Wait(&request,&status);
    MPI_File_close(&matrix1);

    MPI_File matrix2;
    MPI_File_open(MPI_COMM_WORLD,argv[2],MPI_MODE_RDONLY,MPI_INFO_NULL,&matrix2);
    MPI_Datatype new_float2;
    MPI_Type_vector(m[1],m[2],n[2],MPI_FLOAT,&new_float2);
    MPI_Type_commit(&new_float2);
    MPI_File_set_view(matrix2,2*sizeof(int)+sizeof(float)*(dim_rank[1]*n[2]*m[1]+((dim_rank[0]+dim_rank[1])%iter[0])*m[2]),\
                    MPI_FLOAT,new_float2,"native",MPI_INFO_NULL);
    MPI_File_iread(matrix2,&matrix_loc2[0],m[1]*m[2],MPI_FLOAT,&request);
    MPI_Wait(&request,&status);
    MPI_File_close(&matrix2);
  }

  unique_ptr<unique_ptr<float[]>[]> matrix_loc3;      //matrix C for local processor
  double start,end;

  MPI_Barrier(MPI_COMM_WORLD);
  if(iter[1]>1){
    const int iteration[2]={(int)(iter[0]/iter[1]),iter[0]%iter[1]};
    matrix_loc3=std::make_unique<unique_ptr<float[]>[]>(iteration[0]+(iteration[1]>0 ? 1:0));
    for(int i=0;i<iteration[0]+(iteration[1]>0 ? 1:0);i++)
      matrix_loc3[i]=std::make_unique<float[]>(m[0]*m[2]);
    start=MPI_Wtime();
    const int source2=((dim_rank[0]+iter[0]-1+(dim_rank[1]==(iter[1]-1) ?\
                     iter[1]:0))%iter[0])*iter[1]+(dim_rank[1]+1)%iter[1];
    const int dest2=((dim_rank[0]+1+((dim_rank[1]==0) ? (iter[0]-iter[1]):0))\
                    %iter[0])*iter[1]+(dim_rank[1]+iter[1]-1)%iter[1];
    const int source1=(dim_rank[1]+1)%iter[1]+dim_rank[0]*iter[1];
    const int dest1=(dim_rank[1]+iter[1]-1)%iter[1]+dim_rank[0]*iter[1];
    const int source3=((dim_rank[0]+iter[1])%iter[0])*iter[1]+dim_rank[1];
    const int dest3=((dim_rank[0]+iter[0]-iter[1])%iter[0])*iter[1]+dim_rank[1];
    for(int i=0;i<(iteration[0]-1);i+=2){
      for(int j=0;j<(iter[1]-1);j++){
        thread first(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[i][0],m[2],temp_m0,m[1],m[2],iter[2]);
        thread second(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[i][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first.join();
        second.join();
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest2,0,source2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest1,1,source1,1,MPI_COMM_WORLD,&status);
        }
      }
      thread first1(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[i][0],m[2],temp_m0,m[1],m[2],iter[2]);
      thread second1(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[i][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
      first1.join();
      second1.join();
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      }
      for(int j=0;j<(iter[1]-1);j++){
        thread first2(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[i+1][0],m[2],temp_m0,m[1],m[2],iter[2]);
        thread second2(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[i+1][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first2.join();
        second2.join();
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,source2,0,dest2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,source1,1,dest1,1,MPI_COMM_WORLD,&status);
        }
      }
      thread first3(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[i+1][0],m[2],temp_m0,m[1],m[2],iter[2]);
      thread second3(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[i+1][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
      first3.join();
      second3.join();
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      }
    }
    if(iteration[0]%2==1){
      for(int j=0;j<iter[1]-1;j++){
        thread first4(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[iteration[0]-1][0],m[2],temp_m0,m[1],m[2],iter[2]);
        thread second4(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[iteration[0]-1][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first4.join();
        second4.join();
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest2,0,source2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest1,1,source1,1,MPI_COMM_WORLD,&status);
        }
      }
      thread first5(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                  &matrix_loc3[iteration[0]-1][0],m[2],temp_m0,m[1],m[2],iter[2]);
      thread second5(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
            &matrix_loc3[iteration[0]-1][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first5.join();
        second5.join();
    }
    if(iteration[1]>0){
      int k=0;
      while(dim_rank[1]>=(k+1)*iteration[1]){k++;}
      const int source6=((dim_rank[0]+iter[1]-k*iteration[1])%iter[0])*iter[1]+dim_rank[1];
      const int dest6=((dim_rank[0]+iter[0]-iter[1]+k*iteration[1])%iter[0])*iter[1]+dim_rank[1];
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest6,2,source6,2,MPI_COMM_WORLD,&status);
      }
      const int source4=dim_rank[1]<iteration[1]? (dim_rank[1]!=(iteration[1]-1) ?source1:\
                      (rank-(iteration[1]-1))):rank;
      const int dest4=dim_rank[1]<iteration[1]? (dim_rank[1]!=0 ? dest1:\
                    (rank+(iteration[1]-1))):rank;
      const int source5=dim_rank[1]<iteration[1]? (dim_rank[1]!=(iteration[1]-1) ?source2:\
                      (((dim_rank[0]+iteration[1]-1)%iter[0])*iter[1])):\
                      (((dim_rank[0]+1)%iter[0])*iter[1]+dim_rank[1]);
      const int dest5=dim_rank[1]<iteration[1]? (dim_rank[1]!=0 ? dest2:\
                      (iteration[1]-1+((dim_rank[0]+iter[0]-iteration[1]+1)%iter[0])*iter[1])):\
                      (((dim_rank[0]+iter[0]-1)%iter[0])*iter[1]+dim_rank[1]);
      for(int i=0;i<(iteration[1]-1);i++){
        thread first6(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[iteration[0]][0],m[2],temp_m0,m[1],m[2],iter[2]);
        thread second6(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[iteration[0]][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first6.join();
        second6.join();
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest5,0,source5,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest4,1,source4,1,MPI_COMM_WORLD,&status);
        }
        const int color=dim_rank[1]<iteration[1] ? (dim_rank[0]*iter[1]+(dim_rank[1]-i+iteration[1])%iteration[1]):\
                  (dim_rank[0]*iter[1]+dim_rank[1]%iteration[1]);
        int row_rank,row_size;
        MPI_Comm row_comm;
        MPI_Comm_split(MPI_COMM_WORLD,color,rank,&row_comm);
        MPI_Comm_rank(row_comm,&row_rank);
        MPI_Comm_size(row_comm,&row_size);
        const int temp_size=(int)(m[0]/row_size)+(m[0]%row_size>row_rank ? 1:0);
        unique_ptr<float[]>temp_matrix=std::make_unique<float[]>(temp_size*m[2]);
        for(int j=0;j<row_size-1;j++){
          const int send_rank=(row_size+row_rank-j-1)%row_size;
          MPI_Sendrecv(&matrix_loc3[iteration[0]][(send_rank*(int)(m[0]/(row_size))+\
                  (m[0]%row_size>send_rank?send_rank:(m[0]%row_size)))*m[2]],\
                  ((int)(m[0]/row_size)+(m[0]%row_size>send_rank?1:0))*m[2],\
                  MPI_FLOAT,send_rank,0,&temp_matrix[0],temp_size*m[2],MPI_FLOAT,(row_rank+j+1)%row_size,0,row_comm,&status);
          atomic_add(&temp_matrix[0],m[2],&matrix_loc3[iteration[0]][(row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size)))*m[2]],m[2],temp_size,m[2]);
          MPI_Barrier(row_comm);
        }
        memcpy(&temp_matrix[0],&matrix_loc3[(iteration[0]*m[0]+row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size)))*m[2]],sizeof(float)*temp_size*m[2]);
        int* recvcnts=NULL;
        int* disp=NULL;
        if(row_rank==0){
          recvcnts=new int[row_size];
          disp=new int[row_size];
          for(int j=0;j<row_size;j++){
            recvcnts[j]=((int)(m[0]/row_size)+((m[0]%row_size)>j?1:0))*m[2];
            disp[j]=(j*(int)(m[0]/row_size)+((m[0]%row_size)>j?j:(m[0]%row_size)))*m[2];
          }
        }
        MPI_Barrier(row_comm);
        MPI_Igatherv(&temp_matrix[0],temp_size*m[2],MPI_FLOAT,\
                  &matrix_loc3[iteration[0]*m[0]*m[2]],&recvcnts[0],&disp[0],MPI_FLOAT,0,row_comm,&request);
        MPI_Wait(&request,&status);
        if(row_rank!=0){
            memset(&matrix_loc3[iteration[0]*m[0]*m[2]],0,sizeof(float)*m[0]*m[2]);
          }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest4,1,source4,1,MPI_COMM_WORLD,&status);
        }
        delete[] recvcnts,disp;
        MPI_Comm_free(&row_comm);
        }
        thread first7(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[iteration[0]][0],m[2],temp_m0,m[1],m[2],iter[2]);
        thread second7(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[iteration[0]][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
        first7.join();
        second7.join();
        const int color=dim_rank[1]<iteration[1] ? (dim_rank[0]*iter[1]+(dim_rank[1]+1)%iteration[1]):\
                  (dim_rank[0]*iter[1]+dim_rank[1]%iteration[1]);
        int row_rank,row_size;
        MPI_Comm row_comm;
        MPI_Comm_split(MPI_COMM_WORLD,color,rank,&row_comm);
        MPI_Comm_rank(row_comm,&row_rank);
        MPI_Comm_size(row_comm,&row_size);
        const int temp_size=(int)(m[0]/row_size)+(m[0]%row_size>row_rank ? 1:0);
        unique_ptr<float[]>temp_matrix=std::make_unique<float[]>(temp_size*m[2]);
        for(int j=0;j<row_size-1;j++){
          MPI_Sendrecv(&matrix_loc3[iteration[0]][(((row_size+row_rank-j-1)%row_size)*(int)(m[0]/row_size)+\
                  (m[0]%row_size>((row_size+row_rank-j-1)%row_size)?((row_size+row_rank-j-1)%row_size):(m[0]%row_size)))*m[2]],\
                  ((int)(m[0]/row_size)+(m[0]%row_size>((row_size+row_rank-j-1)%row_size)?1:0))*m[2],\
                    MPI_FLOAT,(row_size+row_rank-j-1)%row_size,0,&temp_matrix[0],temp_size*m[2],MPI_FLOAT,(row_rank+j+1)%row_size,0,row_comm,&status);
          atomic_add(&temp_matrix[0],m[2],&matrix_loc3[iteration[0]][(row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size)))*m[2]],m[2],temp_size,m[2]);
          MPI_Barrier(row_comm);
        }
        memcpy(&temp_matrix[0],&matrix_loc3[(iteration[0]*m[0]+(row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size))))*m[2]],sizeof(float)*temp_size*m[2]);
        int* recvcnts=NULL;
        int* disp=NULL;
        if(row_rank==0){
          recvcnts=new int[row_size];
          disp=new int[row_size];
          for(int j=0;j<row_size;j++){
            recvcnts[j]=((int)(m[0]/row_size)+((m[0]%row_size)>j?1:0))*m[2];
            disp[j]=(j*(int)(m[0]/row_size)+((m[0]%row_size)>j?j:(m[0]%row_size)))*m[2];
          }
        }
        MPI_Gatherv(&temp_matrix[0],temp_size*m[2],MPI_FLOAT,\
                  &matrix_loc3[iteration[0]*m[0]*m[2]],recvcnts,disp,MPI_FLOAT,0,row_comm);
        delete[] recvcnts,disp;
        MPI_Comm_free(&row_comm);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    end=MPI_Wtime();
    /*if(temp[0]==0 && temp[1]==0 && temp[2]==0){
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      MPI_Datatype new_float3;
      MPI_Type_vector(m[0],m[2],n[2],MPI_FLOAT,&new_float3);
      MPI_Type_commit(&new_float3);
      for(int i=0;i<iteration[0]+(dim_rank[1]<iteration[1]?1:0);i++){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(((dim_rank[0]+dim_rank[1]+i*iter[1])%iter[0])\
                          *m[2]+dim_rank[0]*n[2]*m[0]),MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_iwrite(matrix3,&matrix_loc3[i][0],m[0]*m[2],MPI_FLOAT,&request);
        MPI_Wait(&request,&status);
      }
      MPI_File_close(&matrix3);
    }
    else{
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      for(int i=0;i<(iteration[0]+(dim_rank[1]<iteration[1]?1:0));i++){
        const int shift=(dim_rank[0]+dim_rank[1]+i*iter[1])%iter[0];
        const int col=m1[2]+(gap[2] > shift ? 1:0);
        const int offset=2*sizeof(int)+sizeof(float)*((gap[2]>shift?(shift*(m1[2]+1)):(gap[2]+shift*m1[2]))+\
                  (gap[0]>dim_rank[0]? (dim_rank[0]*(m1[0]+1)*n[2]):((gap[0]+dim_rank[0]*m1[0])*n[2])));
        for(int j=0;j<row;j++){
          MPI_File_iwrite_at(matrix3,offset+sizeof(float)*j*n[2],&matrix_loc3[i][j*m[2]],col,MPI_FLOAT,&request);
          MPI_Wait(&request,&status);
        }
      }
      MPI_File_close(&matrix3);
    }*/
  }
  else
  {
    matrix_loc3=std::make_unique<unique_ptr<float[]>[]>(NUM_processors);
    for(int i=0;i<NUM_processors;i++)
      matrix_loc3[i]=std::make_unique<float[]>(m[0]*m[2]);
    start=MPI_Wtime();
    const int source=(dim_rank[0]+1)%NUM_processors;
    const int dest=(dim_rank[0]+NUM_processors-1)%NUM_processors;
    for(int i=0;i<(NUM_processors-1);i++){
      thread first8(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[i][0],m[2],temp_m0,m[1],m[2],iter[2]);
      thread second8(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[i][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
      first8.join();
      second8.join();
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest,2,source,2,MPI_COMM_WORLD,&status);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    thread first9(strassen<float>,&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],\
                          &matrix_loc3[NUM_processors-1][0],m[2],temp_m0,m[1],m[2],iter[2]);
    thread second9(strassen<float>,&matrix_loc1[temp_m0*m[1]],m[1],&matrix_loc2[0],m[2],\
              &matrix_loc3[NUM_processors-1][temp_m0*m[2]],m[2],temp_m0,m[1],m[2],iter[2]);
    first9.join();
    second9.join();
    MPI_Barrier(MPI_COMM_WORLD);
    end=MPI_Wtime();
    /*if(temp[0]==0 && temp[1]==0 && temp[2]==0){
      MPI_File matrix3;
      MPI_File_open(MPI_COMM_WORLD,argv[3],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix3);
      MPI_Datatype new_float3;
      MPI_Type_vector(m[0],m[2],n[2],MPI_FLOAT,&new_float3);
      MPI_Type_commit(&new_float3);
      for(int i=0;i<NUM_processors;i++){
        MPI_File_set_view(matrix3,2*sizeof(int)+sizeof(float)*(dim_rank[0]*n[2]*m[0]+((i+dim_rank[0])%NUM_processors)*m[2]),\
                          MPI_FLOAT,new_float3,"native",MPI_INFO_NULL);
        MPI_File_iwrite(matrix3,&matrix_loc3[i][0],m[0]*m[2],MPI_FLOAT,&request);
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
        int col=m1[2]+(gap[2] > k ? 1:0);
        for(int j=0;j<row;j++){
          MPI_File_iwrite_at(matrix3,2*sizeof(int)+sizeof(float)*(offset_loop+j*n[2]),\
                            &matrix_loc3[i][j*m[2]],col,MPI_FLOAT,&request);
          MPI_Wait(&request,&status);
      }
    }
    MPI_File_close(&matrix3);
    }*/
  }
  if(rank==0){
    double time=end-start;
    std::fstream out;
    out.open(argv[3],std::ios::out|std::ios_base::app);
    out<<NUM_processors<<','<<time<<'\n';
    //std::cout<<"Time taken: "<<time<<'\n';
  }
  MPI_Finalize();
  return 0;
}
