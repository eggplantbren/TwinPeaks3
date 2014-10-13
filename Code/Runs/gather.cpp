#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>
#include <algorithm>

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

	// Normalising constant for prior weights
	double logC = -1E300;

	// Prior weights (relative) and 'likelihoods'
	vector<double> logw, logL;

	bool save = false;

	vector<string> directories;
	directories.push_back(string("mcmc_2000/8.1"));
	directories.push_back(string("mcmc_2000/8.2"));
	directories.push_back(string("mcmc_2000/8.3"));
	directories.push_back(string("mcmc_2000/8.4"));
	directories.push_back(string("mcmc_2000/8.5"));
	directories.push_back(string("mcmc_2000/8.6"));
	directories.push_back(string("mcmc_2000/8.7"));

	fstream fout1, fout2;
	if(save)
	{
		fout1.open("logw.txt", ios::out);
		fout2.open("scalars.txt", ios::out);
	}

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
			if(save)
			{
				fout1<<temp1<<endl;
				fout2<<temp2<<' '<<temp3<<endl;
			}

			logC = logsumexp(logC, temp1);
			logw.push_back(temp1);
			logL.push_back(temp2/T1 + temp3/T2);

			k++;
		}

		fin1.close();
		fin2.close();

		cout<<"# Loaded "<<k<<" points."<<endl;
	}

	if(save)
	{
		fout1.close();
		fout2.close();
	}

	// Calculate log(Z)
	double logZ = -1E300;
	for(size_t i=0; i<logL.size(); i++)
		logZ = logsumexp(logZ, (logw[i] - logC) + logL[i]);

	// Calculate partial log(Z)s from each run
	double logw_max = *max_element(logw.begin(), logw.end());
	vector<double> logZ_partial, logC_partial;
	for(size_t i=0; i<logL.size(); i++)
	{
		if(logw[i] == logw_max)
		{
			logZ_partial.push_back(-1E300);
			logC_partial.push_back(-1E300);
		}
		logZ_partial[logZ_partial.size() - 1] = logsumexp(logZ_partial[logZ_partial.size() - 1], logw[i] + logL[i]);
		logC_partial[logC_partial.size() - 1] = logsumexp(logC_partial[logC_partial.size() - 1], logw[i]);
	}
	// Normalise within each run
	for(size_t i=0; i<logZ_partial.size(); i++)
		logZ_partial[i] -= logC_partial[i];
	double max = *max_element(logZ_partial.begin(), logZ_partial.end());
	vector<double> Z(logZ_partial.size());
	double tot = 0.;
	for(size_t i=0; i<logZ_partial.size(); i++)
	{
		Z[i] = exp(logZ_partial[i] - max);
		tot += Z[i];
	}
	for(size_t i=0; i<logZ_partial.size(); i++)
		Z[i] /= tot;
	// Effective number of runs
	double ENR = 0.;
	for(size_t i=0; i<Z.size(); i++)
		ENR += -Z[i]*log(Z[i] + 1E-300);
	ENR = exp(ENR);
	cout<<"# Number of runs = "<<Z.size()<<"."<<endl;
	cout<<"# Effective number of contributing runs = "<<ENR<<"."<<endl;

	
	// Calculate H
	double H = 0.;
	double logP;
	for(size_t i=0; i<logL.size(); i++)
	{
		// Posterior weight
		logP = logw[i] - logC + logL[i] - logZ;
		H += exp(logP)*logL[i];
	}
	H -= logZ;

	// Output
	cout<<setprecision(8);
	cout<<"# ln(Z) = "<<logZ<<"."<<endl;
	cout<<"# H = "<<H<<" nats."<<endl;

	return 0;
}

