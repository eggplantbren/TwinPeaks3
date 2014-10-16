#include "SimpleExample.h"
#include <cassert>
#include <cmath>
#include <RandomNumberGenerator.h>
#include <Utils.h>

using namespace std;
using namespace DNest3;

SimpleExample::SimpleExample()
:Model(2)
,params(200)
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
	int reps = 1 + ((randomU() <= 0.5)?(0):(1 + randInt(99)));
	vector<bool> change(params.size(), false);
	for(int i=0; i<reps; i++)
		change[randInt(params.size())] = true;

	for(int i=0; i<params.size(); i++)
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
	scalars[0] = 0.;
	scalars[1] = 0.;
	for(size_t i=0; i<params.size(); i++)
	{
		scalars[0] += -pow(params[i] - 0.5, 2);
		scalars[1] += -pow(sin(4.*M_PI*params[i]), 2);
	}
}

ostream& operator << (ostream& out, const SimpleExample& m)
{
	for(size_t i=0; i<m.params.size(); i++)
		out<<m.params[i]<<' ';
	return out;
}

