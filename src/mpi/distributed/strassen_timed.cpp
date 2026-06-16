#include <iostream>
#include <fstream>
#include <cmath>
#include <memory>
#include <cstring>
#include <mpi.h>
#include "strassen_mpi.hpp"

using std::unique_ptr;
using strassen_mpi::add;
using strassen_mpi::atomic_add;
using strassen_mpi::atomic_subtract;
using strassen_mpi::C_adjust1;
using strassen_mpi::C_adjust2;
using strassen_mpi::multiply;
using strassen_mpi::nearest_ideal;
using strassen_mpi::strassen;
using strassen_mpi::subtract;

int main(int argc,char** argv){
  MPI_Init(&argc,&argv);
  int rank,NUM_processors;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&NUM_processors);
  int iter[3]={NUM_processors,1,1};
  int temp[3]={0,0,0};
  int n[3]; //the three dimensions whole file and per processor respectively
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
    int threshold=(int)std::pow(2,(int)(-2-std::log2((double)(8*iter[1])/(double)n[1]+(double)(5*iter[0])/(double)(n[0]+n[2]))));
    if(threshold<1){
      threshold=1;
    }
    //std::cout<<"threshold: "<<threshold<<'\n';
    iter[2]=nearest_ideal(n[0],temp[0],iter[0],(n[0]/(4*iter[0]*threshold)));
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

  unique_ptr<float[]> matrix_loc3;      //matrix C for local processor
  double start,end;

  MPI_Barrier(MPI_COMM_WORLD);
  if(iter[1]>1){
    const int iteration[2]={(int)(iter[0]/iter[1]),iter[0]%iter[1]};
    matrix_loc3=std::make_unique<float[]>((iteration[0]+(iteration[1]>0 ? 1:0))*m[0]*m[2]);
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
      int i1=i*m[0]*m[2];
      for(int j=0;j<(iter[1]-1);j++){
        strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[i1],m[2],m[0],m[1],m[2],iter[2]);
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest2,0,source2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest1,1,source1,1,MPI_COMM_WORLD,&status);
        }
      }
      strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[i1],m[2],m[0],m[1],m[2],iter[2]);
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      }
      for(int j=0;j<(iter[1]-1);j++){
        strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[i1+m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,source2,0,dest2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,source1,1,dest1,1,MPI_COMM_WORLD,&status);
        }
      }
      strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[i1+m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest3,2,source3,2,MPI_COMM_WORLD,&status);
      }
    }
    if(iteration[0]%2==1){
      for(int j=0;j<iter[1]-1;j++){
        strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[(iteration[0]-1)*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest2,0,source2,0,MPI_COMM_WORLD,&status);
        }
        for(int k=0;k<m[0];k++){
          MPI_Sendrecv_replace(&matrix_loc1[k*m[1]],m[1],MPI_FLOAT,dest1,1,source1,1,MPI_COMM_WORLD,&status);
        }
      }
      strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[(iteration[0]-1)*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
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
        strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[iteration[0]*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
        for(int k=0;k<m[1];k++){
          MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest5,0,source5,0,MPI_COMM_WORLD,&status);
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
          MPI_Sendrecv(&matrix_loc3[iteration[0]*m[0]*m[2]+(send_rank*(int)(m[0]/(row_size))+\
                  (m[0]%row_size>send_rank?send_rank:(m[0]%row_size)))*m[2]],\
                  ((int)(m[0]/row_size)+(m[0]%row_size>send_rank?1:0))*m[2],\
                  MPI_FLOAT,send_rank,0,&temp_matrix[0],temp_size*m[2],MPI_FLOAT,(row_rank+j+1)%row_size,0,row_comm,&status);
          atomic_add(&temp_matrix[0],m[2],&matrix_loc3[iteration[0]*m[0]*m[2]+(row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size)))*m[2]],m[2],temp_size,m[2]);
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
        strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[iteration[0]*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
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
          MPI_Sendrecv(&matrix_loc3[iteration[0]*m[0]*m[2]+(((row_size+row_rank-j-1)%row_size)*(int)(m[0]/row_size)+\
                  (m[0]%row_size>((row_size+row_rank-j-1)%row_size)?((row_size+row_rank-j-1)%row_size):(m[0]%row_size)))*m[2]],\
                  ((int)(m[0]/row_size)+(m[0]%row_size>((row_size+row_rank-j-1)%row_size)?1:0))*m[2],\
                    MPI_FLOAT,(row_size+row_rank-j-1)%row_size,0,&temp_matrix[0],temp_size*m[2],MPI_FLOAT,(row_rank+j+1)%row_size,0,row_comm,&status);
          atomic_add(&temp_matrix[0],m[2],&matrix_loc3[iteration[0]*m[0]*m[2]+(row_rank*(int)(m[0]/row_size)+(m[0]%row_size>row_rank?row_rank:\
                    (m[0]%row_size)))*m[2]],m[2],temp_size,m[2]);
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
        MPI_File_iwrite(matrix3,&matrix_loc3[i*m[0]*m[2]],m[0]*m[2],MPI_FLOAT,&request);
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
          MPI_File_iwrite_at(matrix3,offset+sizeof(float)*j*n[2],&matrix_loc3[(i*m[0]+j)*m[2]],col,MPI_FLOAT,&request);
          MPI_Wait(&request,&status);
        }
      }
      MPI_File_close(&matrix3);
    }*/
  }
  else
  {
    matrix_loc3=std::make_unique<float[]>(NUM_processors*m[0]*m[2]);
    start=MPI_Wtime();
    const int source=(dim_rank[0]+1)%NUM_processors;
    const int dest=(dim_rank[0]+NUM_processors-1)%NUM_processors;
    for(int i=0;i<(NUM_processors-1);i++){
      strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[i*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
      for(int k=0;k<m[1];k++){
        MPI_Sendrecv_replace(&matrix_loc2[k*m[2]],m[2],MPI_FLOAT,dest,2,source,2,MPI_COMM_WORLD,&status);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    }
    strassen(&matrix_loc1[0],m[1],&matrix_loc2[0],m[2],&matrix_loc3[(NUM_processors-1)*m[0]*m[2]],m[2],m[0],m[1],m[2],iter[2]);
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
        MPI_File_iwrite(matrix3,&matrix_loc3[i*m[0]*m[2]],m[0]*m[2],MPI_FLOAT,&request);
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
                            &matrix_loc3[(i*m[0]+j)*m[2]],col,MPI_FLOAT,&request);
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
