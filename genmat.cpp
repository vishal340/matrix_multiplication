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

	ofstream ofile(output_file,ofstream::out|ofstream::trunc);
	if(!ofile)
		throw runtime_error("Could not open file "+string(output_file));
        ofile << n1 <<"\t"<<n2<< "\n";
	for(int i=0; i<n1; i++)
	{
		for(int j=0; j<(n2-1); j++)
		{
			arr[i*n1+j] = ((float) rand())/(RAND_MAX/10)-5;
			ofile << fixed << setprecision(3) << (float)((float)(arr[i*n1+j]*10))/10 << '\t';
		}
		arr[i*n1+n2-1] = ((float) rand())/(RAND_MAX/10)-5;
		ofile << fixed << setprecision(3) << (float)((float)(arr[i*n1+n2-1]*10))/10 << '\n';
	}

	delete arr;
	ofile.close();
	return 0;
}
