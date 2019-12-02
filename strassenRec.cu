#include<cstdio>
#include<fstream>
#include<cmath>
#include<cuda.h>

int threshold=256;
int xthread=32;



__global__ void  multiply(float* A,float* B,float* C,int jump,int jump1,int jump2,int iter)
{
        __shared__ float A1[32][32],B1[32][32];
	int posy=blockIdx.y*blockDim.y+threadIdx.y;
	int posx=blockIdx.x*blockDim.x+threadIdx.x;
	int row=posy*jump+threadIdx.x;
	int col=posx+threadIdx.y*jump1;
        int place=posy*jump2+posx;
        for(int i=0;i<iter;i++)
        {
                A1[threadIdx.y][threadIdx.x]=A[row];
                B1[threadIdx.y][threadIdx.x]=B[col];
                __syncthreads();
                for(int i=0;i<blockDim.x;i++)
                        C[place]+=A1[threadIdx.y][i]*B1[i][threadIdx.x];
                row+=blockDim.x;
                col+=blockDim.y*jump1;
		__syncthreads();
        }
}

//non_square matrix multiplication

__global__ void kernel1(float* A,int jump1,int n1,float* B,int jump2,int n2,float* C,int jump3,float* D,int jump4,int m1,int m2,int m3)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y+threadIdx.y;
	if(col<m2 && row<m1)
		C[col+jump3*row]=A[col+jump1*row]+A[n1+col+jump1*row];
	if(col<m3 && row<m2)
		D[col+jump4*row]=B[col+jump2*row]+B[n2+col+jump2*row];
}
__global__ void kernel7(float* A,int jump1,int n1,float* B,int jump2,int n2,float* C,int jump3,float* D,int jump4,int m1,int m2,int m3)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y+threadIdx.y;
	if(col<m2 && row<m1)
		C[col+jump3*row]=A[col+jump1*row]-A[n1+col+jump1*row];
	if(col<m3 && row<m2)
		D[col+jump4*row]=B[col+jump2*row]+B[n2+col+jump2*row];
}

__global__ void kernel2(float* A,int jump,int n1,int n2,float* B,int jump1,float* C,int jump2,int n3,int m2,int m3)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	if(col<m3)
	{
	if(3*blockIdx.y<gridDim.y)
	{
		int point=col+(blockIdx.y * blockDim.y + threadIdx.y)*jump;
		A[point+n2]+=A[point];
	}
	else if(3*blockIdx.y<2*gridDim.y)
	{
		int point=col+((blockIdx.y-gridDim.y/3) * blockDim.y + threadIdx.y)*jump;
		A[point]+=A[point+n1];
		A[point+n1]=0;
	}
	}
	if(col<m2 && 3*blockIdx.y>=2*gridDim.y)
	{
		int row=(blockIdx.y-2*gridDim.y/3) * blockDim.y + threadIdx.y;
		B[col+jump1*row]=C[col+jump2*row]+C[col+jump2*row+n3];
	}
}

__global__ void kernel3(float* A,int jump,int n,float* B,int jump1,float* C,int jump2,int n1,int m2,int m3)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y + threadIdx.y;
	if(col<m2)
		A[col+jump*row+n]-=A[col+jump*row];
	if(col<m3)
		B[col+jump1*row]=C[col+jump2*row]+C[col+jump2*row+n1];
}

__global__ void kernel4(float* A,int jump,int n1,float* B,int jump1,float* C,int jump2,float* D,int jump3,int n2,int m1,int m2)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x,row;
	if(2*blockIdx.y<gridDim.y)
	{
		row=blockIdx.y * blockDim.y + threadIdx.y;
		if(row<m1)
		{
			A[col+jump*row]+=B[col+jump1*row];
			A[col+jump*row+n1]+=B[col+jump1*row];
			B[col+jump1*row]=0;
		}
	}
	else
	{
		row=(blockIdx.y-gridDim.y/2) * blockDim.y + threadIdx.y;
		if(row<m2)
			C[col+jump2*row]=D[col+jump3*row+n2]-D[col+jump3*row];
	}
}

__global__ void kernel5(float* A,int jump,int n1,float* B,int jump1)
{
	int point=(blockIdx.x * blockDim.x+threadIdx.x)+(blockIdx.y * blockDim.y + threadIdx.y)*jump;
	int point1=(blockIdx.x * blockDim.x+threadIdx.x)+(blockIdx.y * blockDim.y + threadIdx.y)*jump1;
	A[point]+=B[point1];
	A[point+n1]+=B[point1];
}

__global__ void kernel6(float* A,int jump,float* B,int jump1,int n1,float* C,int jump2,float* D,int jump3,int n2,int m1,int m2,int m3)
{
	int col=blockIdx.x*blockDim.x+threadIdx.x;
	if(col<m2 && 3*blockIdx.y<gridDim.y)
	{
		int row=blockIdx.y*blockDim.y+threadIdx.y;
		if(row<m1)
		A[col+row*jump]=0;
	}
	if(col<m3 && 3*blockIdx.y>=gridDim.y && 3*blockIdx.y<2*gridDim.y)
	{
		int row=(blockIdx.y-gridDim.y/3) * blockDim.y + threadIdx.y;
		if(row<m1)
		B[col+row*jump1]-=B[col+row*jump1+n1];
	}
	if(col<m3 && 3*blockIdx.y>2*gridDim.y)
	{
		int row=(blockIdx.y-(2*gridDim.y/3)) * blockDim.y + threadIdx.y;
		if(row<m2)
		C[col+jump2*row]=D[col+jump3*row]-D[col+jump3*row+n2];
	}
}

void strassen(float* A,int jump,float* B,int jump1,float* C,int jump2,int n1,int n2,int n3,float* temp1,float* temp2,int n,int n_)
{
  if(n1<=threshold || n2<=threshold || n3<=threshold)
  {
    multiply <<<dim3(n3/xthread,n1/xthread),dim3(xthread,xthread)>>> (A,B,C,jump,jump1,jump2,n2/xthread);
  }
  else
  {
    n1/=2;n2/=2;n3/=2;
    n/=2;n_/=2;
    //M1
    kernel1 <<<dim3(n_/xthread,n/xthread),dim3(xthread,xthread)>>> (A,jump,jump*n1+n2,B,jump1,jump1*n2+n3,temp1,n2,temp2,n3,n1,n2,n3);   //temp1=A11+A22 //temp2=B11+B22
    strassen(temp1,n2,temp2,n3,C,jump2,n1,n2,n3,temp1+n1*n2,temp2+n2*n3,n,n_);   //C11=temp1*temp2
    //M6
    kernel7 <<<dim3(n_/xthread,n/xthread),dim3(xthread,xthread)>>> (A+jump*n1,jump,-jump*n1,B,jump1,n3,temp1,n2,temp2,n3,n1,n2,n3);   //temp1=A21-A11//temp2=B11+B12
    strassen(temp1,n2,temp2,n3,C+jump2*n1+n3,jump2,n1,n2,n3,temp1+n1*n2,temp2+n2*n3,n,n_);   //C22=temp1*temp2
    //M7
    kernel7 <<<dim3(n_/xthread,n/xthread),dim3(xthread,xthread)>>> (A+n2,jump,jump*n1,B+jump1*n2,jump1,n3,temp1,n2,temp2,n3,n1,n2,n3);//temp1=A12-A22 //temp2=B21+B22
    strassen(temp1,n2,temp2,n3,C+n3,jump2,n1,n2,n3,temp1+n1*n2,temp2+n2*n3,n,n_);  //C12+=temp1*temp2
    kernel2 <<<dim3(n_/xthread,(3*n1)/xthread),dim3(xthread,xthread)>>> (C,jump2,n3,jump2*n1+n3,temp1,n2,A+jump*n1,jump,n2,n2,n3);  //C22+=C11//C11+=C12 //temp1=A21+A22 //C12=0
    //M2
    strassen(temp1,n2,B,jump1,C+jump2*n1,jump2,n1,n2,n3,temp1+n1*n2,temp2+n2*n3,n,n_);  //C21=temp1*B11
    kernel3 <<<dim3(n_/xthread,n1/xthread),dim3(xthread,xthread)>>> (C+jump2*n1,jump2,n3,temp1,n2,A,jump,n2,n2,n3);       //C22-=C21 //temp1=A11+A12
    //M5
    strassen(temp1,n2,B+jump1*n2+n3,jump1,C+n3,jump2,n1,n2,n3,temp1+n1*n2,temp2+n2*n3,n,n_);  //C12=temp1*B22
    kernel6<<<dim3(n_/xthread,(3*n)/xthread),dim3(xthread,xthread)>>>(temp1,n2,C,jump2,n3,temp2,n3,B+n3,jump1,jump1*n2,n1,n2,n3); //C11-=C12 //temp2=B12-B22 //temp1=0
    //M3
    strassen(A,jump,temp2,n3,temp1,n3,n1,n2,n3,temp1+n1*n3,temp2+n2*n3,n,n_);  //temp1=A11*temp2
    kernel4 <<<dim3(n3/xthread,(2*n)/xthread),dim3(xthread,xthread)>>> (C+n3,jump2,jump2*n1,temp1,n2,temp2,n3,B,jump1,jump1*n2,n1,n2);  //C12=C12+temp1//C22=C22+temp1 //temp2=B21-B11 //temp1=0
    //M4
    strassen(A+jump*n1+n2,jump,temp2,n3,temp1,n3,n1,n2,n3,temp1+n1*n3,temp2+n2*n3,n,n_);  //temp1=A22*temp2
    kernel5 <<<dim3(n3/xthread,n1/xthread),dim3(xthread,xthread)>>> (C,jump2,jump2*n1,temp1,n2);                          //C21=C21+temp1//C11=C11+temp1
  }
}

//Square matrix multiplication

__global__ void kernel1_sq(float* A,int jump1,int n1,float* B,int jump2,int n2,float* C,float* D,int jump) 
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y+threadIdx.y;
	C[col+jump*row]=A[col+jump1*row]+A[n1+col+jump1*row];
	D[col+jump*row]=B[col+jump2*row]+B[n2+col+jump2*row];
}
__global__ void kernel7_sq(float* A,int jump1,int n1,float* B,int jump2,int n2,float* C,float* D,int jump)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y+threadIdx.y;
	C[col+jump*row]=A[col+jump1*row]-A[n1+col+jump1*row];
	D[col+jump*row]=B[col+jump2*row]+B[n2+col+jump2*row];
}

__global__ void kernel2_sq(float* A,int jump,int n1,int n2,float* B,int jump1,float* C,int jump2,int n3)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	if(3*blockIdx.y<gridDim.y)
	{
		int point=col+(blockIdx.y * blockDim.y + threadIdx.y)*jump;
		A[point+n2]+=A[point];
	}
	else if(3*blockIdx.y<2*gridDim.y)
	{
		int point=col+((blockIdx.y-gridDim.y/3) * blockDim.y + threadIdx.y)*jump;
		A[point]+=A[point+n1];
		A[point+n1]=0;
	}
	else
	{
		int row=(blockIdx.y-2*gridDim.y/3) * blockDim.y + threadIdx.y;
		B[col+jump1*row]=C[col+jump2*row]+C[col+jump2*row+n3];
	}
}

__global__ void kernel3_sq(float* A,int jump,int n,float* B,int jump1,float* C,int jump2,int n1)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x;
	int row=blockIdx.y * blockDim.y + threadIdx.y;
	A[col+jump*row+n]-=A[col+jump*row];
	B[col+jump1*row]=C[col+jump2*row]+C[col+jump2*row+n1];
}

__global__ void kernel4_sq(float* A,int jump,int n1,float* B,int jump1,float* C,int jump2,float* D,int jump3,int n2)
{
	int col=blockIdx.x * blockDim.x+threadIdx.x,row;
	if(2*blockIdx.y<gridDim.y)
	{
		row=blockIdx.y * blockDim.y + threadIdx.y;
		A[col+jump*row]+=B[col+jump1*row];
		A[col+jump*row+n1]+=B[col+jump1*row];
		B[col+jump1*row]=0;
	}
	else
	{
		row=(blockIdx.y-gridDim.y/2) * blockDim.y + threadIdx.y;
		C[col+jump2*row]=D[col+jump3*row+n2]-D[col+jump3*row];
	}
}

__global__ void kernel6_sq(float* A,int jump,float* B,int jump1,int n1,float* C,int jump2,float* D,int jump3,int n2)
{
	int col=blockIdx.x*blockDim.x+threadIdx.x;
	if(3*blockIdx.y<gridDim.y)
	{
		int row=blockIdx.y*blockDim.y+threadIdx.y;
		A[col+row*jump]=0;
	}
	else if(3*blockIdx.y<2*gridDim.y)
	{
		int row=(blockIdx.y-gridDim.y/3) * blockDim.y + threadIdx.y;
		B[col+row*jump1]-=B[col+row*jump1+n1];
	}
	else
	{
		int row=(blockIdx.y-(2*gridDim.y/3)) * blockDim.y + threadIdx.y;
		C[col+jump2*row]=D[col+jump3*row]-D[col+jump3*row+n2];
	}
}
void strassen_sq(float* A,int jump,float* B,int jump1,float* C,int jump2,int n,float* temp1,float* temp2,int block_len)
{
  if(n<=threshold)
	multiply<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(A,B,C,jump,jump1,jump2,n/xthread);
  else
  {
	n/=2;
	block_len/=2;
     	kernel1_sq<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(A,jump,jump*n+n,B,jump1,jump1*n+n,temp1,temp2,n);
	strassen_sq(temp1,n,temp2,n,C,jump2,n,temp1+n*n,temp2+n*n,block_len);
        kernel7_sq<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(A+jump*n,jump,-jump*n,B,jump1,n,temp1,temp2,n);
	strassen_sq(temp1,n,temp2,n,C+jump2*n+n,jump2,n,temp1+n*n,temp2+n*n,block_len);
     	kernel7_sq<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(A+n,jump,jump*n,B+jump1*n,jump1,n,temp1,temp2,n);
	strassen_sq(temp1,n,temp2,n,C+n,jump2,n,temp1+n*n,temp2+n*n,block_len);
        kernel2_sq<<<dim3(block_len,3*block_len),dim3(xthread,xthread)>>>(C,jump2,n,jump2*n+n,temp1,n,A+jump*n,jump,n);
	strassen_sq(temp1,n,B,jump1,C+jump2*n,jump2,n,temp1+n*n,temp2+n*n,block_len);
	kernel3_sq<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(C+jump2*n,jump2,n,temp1,n,A,jump,n);
	strassen_sq(temp1,n,B+jump1*n+n,jump1,C+n,jump2,n,temp1+n*n,temp2+n*n,block_len);
        kernel6_sq<<<dim3(block_len,3*block_len),dim3(xthread,xthread)>>>(temp1,n,C,jump2,n,temp2,n,B+n,jump1,jump1*n);
	strassen_sq(A,jump,temp2,n,temp1,n,n,temp1+n*n,temp2+n*n,block_len);
 	kernel4_sq<<<dim3(block_len,2*block_len),dim3(xthread,xthread)>>>(C+n,jump2,jump2*n,temp1,n,temp2,n,B,jump1,jump1*n);
	strassen_sq(A+jump*n+n,jump,temp2,n,temp1,n,n,temp1+n*n,temp2+n*n,block_len);
	kernel5<<<dim3(block_len,block_len),dim3(xthread,xthread)>>>(C,jump2,jump2*n,temp1,n);
  }
}
 

int nearest_ideal(int &n,int &temp)
{
  int temp1=(xthread-n%xthread)%xthread;
  int pow=1;
  n+=temp1;
  int m=n/xthread;
  while(m>threshold/xthread)
  {
    if(m%2==1)
    {
        temp+=pow;
        m++;
    }
    m/=2;
    pow*=2;
  }
  n+=temp*xthread;
  temp=temp*xthread;
  temp+=temp1;
  return pow;
}

int main(int argc,char** argv)
{
  std::ifstream in(argv[1]);
  std::ifstream in1(argv[2]);
  std::ofstream out(argv[3]);
  float *A,*B,*C;
  int n1,n2,n3;
  int temp1=0,temp2=0,temp3=0;
  in>>n1>>n2;
  in1>>n2>>n3;
  out<<n1<<'\t'<<n3<<'\n';
  int power=nearest_ideal(n1,temp1);
  power=std::min(power,nearest_ideal(n2,temp2));
  power=std::min(power,nearest_ideal(n3,temp3));
  float factor=0;
  for(int i=power;i>1;i/=2)
     factor+=1/(float)(i*i);
  A=new float[n1*n2];
  B=new float[n2*n3];
  C=new float[n1*n3];
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
  int size_temp1,size_temp2,n_,n;
  n=n2>n3?n2:n3;
  n_=n1>n2?n1:n2;
  size_temp1=(int)(n1*n_)*factor;
  size_temp2=(int)(n2*n3)*factor;
  float *d_A, *d_B, *d_C,*temp1_,*temp2_;
  cudaMalloc( (void **) &d_A, sizeof(float)*n1*n2);
  cudaMalloc( (void **) &d_B, sizeof(float)*n2*n3);
  cudaMalloc( (void **) &d_C, sizeof(float)*n1*n3);
  cudaMalloc( (void **) &temp1_,sizeof(float)*size_temp1);
  cudaMalloc( (void **) &temp2_,sizeof(float)*size_temp2);
  //copy from host to device
  cudaMemcpy (d_A, A, sizeof(float)*n1*n2, cudaMemcpyHostToDevice);
  cudaMemcpy (d_B, B, sizeof(float)*n2*n3, cudaMemcpyHostToDevice);
  cudaMemset(d_C,0,sizeof(float)*n1*n3);
  cudaMemset(temp1_,0,sizeof(float)*size_temp1);
  cudaMemset(temp2_,0,sizeof(float)*size_temp2);
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&start);
  if(n1!=n2 || n3!=n2)
  	strassen(d_A,n2,d_B,n3,d_C,n3,n1,n2,n3,temp1_,temp2_,n,n_);
  else
  	strassen_sq(d_A,n1,d_B,n1,d_C,n1,n1,temp1_,temp2_,n1/xthread);
  cudaDeviceSynchronize(); 
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&end);

  cudaMemcpy (C, d_C, sizeof(float)*n1*n3, cudaMemcpyDeviceToHost);
  printf("Error %s \n",cudaGetErrorString(cudaGetLastError()));
  double time_taken = (end.tv_nsec-start.tv_nsec)+1e+9*(end.tv_sec-start.tv_sec);
  printf("StrassenRec - Time taken: %f\n",time_taken);
  for(int i=0;i<n1-temp1;i++)
  {
    for(int j=0;j<n3-temp3;j++)
      out<<C[i*n3+j]<<'\t';
    out<<'\n';
  }
  std::ofstream ofile;
  ofile.open(argv[4],std::ios_base::app);
  ofile<<"strassenRec - Time taken (ns): "<<time_taken<<"\n";
  ofile.close();

  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_C);
  cudaFree(temp1_);
  cudaFree(temp2_);
  delete C,A,B;
}
