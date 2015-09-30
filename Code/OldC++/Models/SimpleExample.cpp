#include "SimpleExample.h"
#include <cassert>
#include <cmath>
#include <RandomNumberGenerator.h>
#include <Utils.h>

using namespace std;
using namespace DNest3;

SimpleExample::SimpleExample()
:Model(2)
,params(20)
{

}

SimpleExample::~SimpleExample()
{

}

void SimpleExample::from_prior()
{
	for(size_t i=0; i<params.size(); i++)
		params[i] = randomU();
	compute_scalars();
}

double SimpleExample::perturb()
{
	int reps = 1 + ((randomU() <= 0.5)?(0):(1 + randInt(9)));
	vector<bool> change(params.size(), false);
	for(int i=0; i<reps; i++)
		change[randInt(params.size())] = true;

	for(size_t i=0; i<params.size(); i++)
	{
		if(change[i])
		{
			params[i] += randh();
			wrap(params[i], 0., 1.);
		}
	}
	compute_scalars();

	return 0.;
}

void SimpleExample::compute_scalars()
{
	double u = 0.01;
	double v = 0.1;
	double C = log(1.0/sqrt(2*M_PI));
	double logl1 = 0;
	double logl2 = 0;

	for(size_t i=0; i<params.size(); i++)
	{
		logl1 += C - log(u) - 0.5*pow((params[i] - 0.0)/u, 2);
		logl2 += C - log(v) - 0.5*pow(params[i]/v, 2);
	}
	logl1 += log(100.0);

	scalars[0] = logsumexp(logl1, logl2);
	scalars[1] = 0.;
	for(size_t i=0; i<params.size(); i++)
		scalars[1] += -pow(sin(4.*M_PI*params[i]), 2);
}

ostream& operator << (ostream& out, const SimpleExample& m)
{
	for(size_t i=0; i<m.params.size(); i++)
		out<<m.params[i]<<' ';
	return out;
}

