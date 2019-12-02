#include <iostream>
#include <string.h>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
	int n1,n2;

	char output_file[5+strlen(argv[1])+strlen(argv[2])];
	strcpy(output_file,"in_");
	strcat(output_file, (const char*)argv[1]);
	strcat(output_file,"_");
	strcat(output_file,(const char*)argv[2]);

    istringstream iss(argv[1]);
	istringstream iss1(argv[2]);
    iss>>n1;
	iss1>>n2;
	float *arr = new float[n1*n2]();

	ofstream ofile;
	ofile.open(output_file,ios::out|ios::binary);
	if(!ofile)
		throw runtime_error("Could not open file "+string(output_file));

        ofile.write((char*)&n1,sizeof(int));
		ofile.seekp(sizeof(int),ios::beg);
		ofile.write((char*)&n2,sizeof(int));
	for(int i=0; i<n1; i++)
	{
		for(int j=0; j<n2; j++)
		{
			arr[i*n1+j] = ((float) rand())/(RAND_MAX/10)-5;
		}
	}
	ofile.seekp(sizeof(int)*2,ios::beg);
	ofile.write((char*)arr,sizeof(float)*n1*n2);
	delete arr;
	ofile.close();
	return 0;
}
