#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>

using namespace std;

void multiply(int size, float *matA, float* matB, float *prod)
{
	for(int i=0; i<size; i++)
	{
		for(int k=0; k<size; k++)
		{
			for(int j=0; j<size; j++)
			{
				prod[i*size+j] += matA[i*size+k]*matB[k*size+j];
			}
		}
	}
}

int main(int argc, char *argv[])
{
	/*if(argc != 4 )
	{
		cout<<"Invalid input arguments.\n";
		return 0;
	}*/
	struct timeval start, stop;
	char *ipfile = argv[2];
	char *outfile = argv[3];
	//char *logfile = argv[3];
	int size=atoi(argv[1]);

	ifstream ifile(ipfile,ifstream::in);
	if(!ifile)
	{
		cout<<"Can't open input file.\n";
		return 0;
	}
	float *mat = new float[size*size];
	int n1,n2;
	ifile>>n1>>n2;
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)
			ifile>>mat[i*size+j];
	}
	ifile.close();

	float *prod = new float[size*size]();
	gettimeofday(&start, 0);
	multiply(size, mat, mat, prod);
	gettimeofday(&stop, 0);

	ofstream ofile(outfile,ofstream::out);
	for(int i=0; i<size; i++)
	{
		for(int j=0; j<size; j++)
			ofile<<prod[i*size+j]<<" ";
		ofile<<"\n";
	}
	ofile.close();

        double time_taken = (1000000.0*(stop.tv_sec-start.tv_sec) + stop.tv_usec-start.tv_usec);

	//ofile.open(logfile,ios_base::app);
        //ofile<<"SimpleIKJ - Time taken (us): "<< time_taken<<"\n";
	//ofile.close();
	cout<<"Time taken: "<<(stop.tv_usec-start.tv_usec)+1e+6*(stop.tv_sec-start.tv_sec);
	delete prod;
	delete mat;
	return 0;
}
