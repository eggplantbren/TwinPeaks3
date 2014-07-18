#include "Model.h"
#include <RandomNumberGenerator.h>
#include <Utils.h>

using namespace DNest3;
using namespace std;

Model::Model(int num_scalars)
:num_scalars(num_scalars), scalars(num_scalars), tiebreakers(num_scalars)
{

}

Model::~Model()
{

}

void Model::from_prior_tiebreakers()
{
	for(size_t i=0; i<tiebreakers.size(); i++)
		tiebreakers[i] = randomU();
}

double Model::perturb_tiebreakers()
{
	int which = randInt(tiebreakers.size());
	tiebreakers[which] += randh();
	wrap(tiebreakers[which], 0., 1.);
	return 0.;
}

bool Model::is_above(const vector< vector<double> >& threshold) const
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

