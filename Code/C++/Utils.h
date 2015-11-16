#ifndef _Utils_
#define _Utils_

// Useful functions, copied from DNest3

#include <vector>
#include <algorithm>

double mod(double y, double x);
int mod(int y, int x);
void wrap(double& x, double min, double max);
double logsumexp(double* logv, int n);
double logsumexp(const std::vector<double>& logv);
double logsumexp(double a, double b);
double logdiffexp(double a, double b);

// Argsort from
// http://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> argsort(const std::vector<T>& v)
{
	// initialize original index locations
	std::vector<size_t> idx(v.size());
	for(size_t i=0; i<idx.size(); i++)
		idx[i] = i;

	// sort indexes based on comparing values in v
	std::sort(idx.begin(), idx.end(),
		[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

	return idx;
}

#endif

