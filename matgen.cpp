#include<iostream>
#include<cstdlib>
#include"mpi.h"
using namespace std;

int main(int argc,char** argv)
{
MPI_Init(&argc,&argv);
int rank,NUM_processors;
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
MPI_Comm_size(MPI_COMM_WORLD,&NUM_processors);
int row_len=atoi(argv[3]);
int col_len=atoi(argv[2]);
int stride=row_len/NUM_processors;
MPI_File matrix;
MPI_File_open(MPI_COMM_WORLD,argv[1],MPI_MODE_WRONLY|MPI_MODE_CREATE,MPI_INFO_NULL,&matrix);
MPI_Status status;
if(rank==0){
    MPI_File_write(matrix,&col_len,1,MPI_INT,&status);
    MPI_File_write(matrix,&row_len,1,MPI_INT,&status);
}
MPI_Datatype new_int;
MPI_Type_vector(col_len,stride,row_len,MPI_FLOAT,&new_int);
MPI_Type_commit(&new_int);
MPI_File_set_view(matrix,sizeof(int)*rank*stride,MPI_FLOAT,new_int,"native",MPI_INFO_NULL);

float **y=new float*[col_len];
for(int i=0;i<col_len;i++)
y[i]=new float[stride];
for(int i=0;i<col_len;i++)
{
for(int j=0;j<stride;j++)
y[i][j]=(float)((float) rand())/(RAND_MAX/10)-5;
}
MPI_File_write_all_begin(matrix,y,col_len*stride,MPI_FLOAT);
MPI_File_write_all_end(matrix,y,MPI_STATUS_IGNORE);
MPI_File_close(&matrix);
MPI_Finalize();
return 0;
}
