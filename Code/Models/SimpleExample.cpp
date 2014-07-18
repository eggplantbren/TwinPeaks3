#include "SimpleExample.h"
#include <cassert>
#include <cmath>
#include <RandomNumberGenerator.h>
#include <Utils.h>

using namespace std;
using namespace DNest3;

SimpleExample::SimpleExample()
:params(100)
,scalars(2)
,tiebreakers(2)
{

}

void SimpleExample::from_prior()
{
	for(size_t i=0; i<params.size(); i++)
		params[i] = randomU();
	compute_scalars();
}

void SimpleExample::from_prior_tiebreakers()
{
	for(size_t i=0; i<tiebreakers.size(); i++)
		tiebreakers[i] = randomU();
}

double SimpleExample::perturb_tiebreakers()
{
	int reps = 1 + ((randomU() <= 0.5)?(0):(1 + randInt(tiebreakers.size()-1)));

	for(int i=0; i<reps; i++)
	{
		int which = randInt(tiebreakers.size());
		tiebreakers[which] += randh();
		wrap(tiebreakers[which], 0., 1.);
	}

	return 0.;
}

double SimpleExample::perturb()
{
	int reps = 1 + ((randomU() <= 0.5)?(0):(1 + randInt(99)));

	for(int i=0; i<reps; i++)
	{
		int which = randInt(params.size());
		params[which] += randh();
		wrap(params[which], 0., 1.);
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

bool SimpleExample::is_above(const vector< vector<double> >& threshold) const
{
	assert(scalars.size() == threshold.size());
	for(size_t i=0; i<scalars.size(); i++)
	{
		if((scalars[i] < threshold[i][0]) ||
			((scalars[i] == threshold[i][0] && tiebreakers[i] < threshold[i][1])))
			return false;
	}
	return true;
}

ostream& operator << (ostream& out, const SimpleExample& m)
{
	for(size_t i=0; i<m.params.size(); i++)
		out<<m.params[i]<<' ';
	return out;
}

