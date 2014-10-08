#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main()
{
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
			k++;
		}

		fin1.close();
		fin2.close();

		cout<<"# Processed "<<k<<" points."<<endl;
	}

	fout1.close();
	fout2.close();

	return 0;
}

