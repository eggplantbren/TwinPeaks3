#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;


double logsumexp(double x1, double x2)
{
	if(x1 > x2)
		return x1 + log(1. + exp(x2 - x1));
	return x2 + log(1. + exp(x1 - x2));
}


int main()
{
	// Temperatures
	double T1 = 0.3; double T2 = 0.3;
	double log_numerator = -1E300; double log_denominator = -1E300;

	vector<string> directories;
	directories.push_back(string("10.1"));
	directories.push_back(string("10.2"));
	directories.push_back(string("10.3"));
	directories.push_back(string("10.4"));
	directories.push_back(string("10.5"));
	directories.push_back(string("10.6"));
	directories.push_back(string("10.7"));
	directories.push_back(string("10.8"));

	fstream fout1("logw.txt", ios::out);
	fstream fout2("scalars.txt", ios::out);

	int k = 0;
	for(size_t i=0; i<directories.size(); i++)
	{
		string filename1 = directories[i] + string("/logw.txt");
		fstream fin1(filename1.c_str(), ios::in);

		string filename2 = directories[i] + string("/scalars.txt");
		fstream fin2(filename2.c_str(), ios::in);

		double temp1, temp2, temp3;
		while(fin1>>temp1 && fin2>>temp2 && fin2>>temp3)
		{
			fout1<<temp1<<endl;
			fout2<<temp2<<' '<<temp3<<endl;

			log_numerator = logsumexp(log_numerator, temp1 + temp2/T1 + temp3/T2);
			log_denominator = logsumexp(log_denominator, temp1);
			k++;
		}

		fin1.close();
		fin2.close();

		cout<<"# Processed "<<k<<" points."<<endl;
	}

	fout1.close();
	fout2.close();

	cout<<setprecision(8)<<"# ln(Z) = "<<(log_numerator - log_denominator)<<endl;

	return 0;
}

