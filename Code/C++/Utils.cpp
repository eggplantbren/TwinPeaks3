#include "Utils.h"
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

double mod(double y, double x)
{
	if(x <= 0)
		cerr<<"Warning in mod(double, double) (Utils.cpp)"<<endl;
	return (y/x - floor(y/x))*x;
}

void wrap(double& x, double min, double max)
{
	x = mod(x - min, max - min) + min;
}

int mod(int y, int x)
{
	if(x <= 0)
		cerr<<"Warning in mod(int, int) (Utils.cpp)"<<endl;
	if(y >= 0)
		return y - (y/x)*x;
	else
		return (x-1) - mod(-y-1, x);
}

double logsumexp(double* logv, int n)
{
	if(n<=1)
		cerr<<"Warning in logsumexp"<<endl;
	double max = logv[0];
	for(int i=1; i<n; i++)
		if(logv[i] > max)
			max = logv[i];
	double answer = 0;
	// log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
	//			= max(logf) + log(sum(exp(logf - max(logf)))
	for(int i=0; i<n; i++)
		answer += exp(logv[i] - max);
	answer = max + log(answer);
	return answer;
}

double logsumexp(const vector<double>& logv)
{
	int n = static_cast<int>(logv.size());
	//if(n<=1)
	//	cout<<"Warning in logsumexp"<<endl;
	double max = *max_element(logv.begin(), logv.end());
	double answer = 0;
	// log(sum(exp(logf)) 	= log(sum(exp(logf - max(logf) + max(logf)))
	//			= max(logf) + log(sum(exp(logf - max(logf)))
	for(int i=0; i<n; i++)
		answer += exp(logv[i] - max);
	answer = max + log(answer);
	return answer;
}

double logsumexp(double a, double b)
{
	double x[2] = {a,b};
	return logsumexp(x, 2);
}

double logdiffexp(double a, double b)
{
	if(a <= b)
		cerr<<"# Error in logdiffexp."<<endl;
	double biggest = a;
	a -= biggest;
	b -= biggest;
	return log(exp(a) - exp(b)) + biggest;
}

