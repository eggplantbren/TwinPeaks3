#include "MyModel.h"
#include <cmath>
#include <RandomNumberGenerator.h>
#include <Utils.h>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:params(100)
,scalars(2)
{

}

void MyModel::from_prior()
{
	for(size_t i=0; i<params.size(); i++)
		params[i] = randomU();
	compute_scalars();
}

double MyModel::perturb()
{
	int reps = 1 + (randomU() <= 0.5)?(0):(1 + randInt(9));

	for(int i=0; i<reps; i++)
	{
		int which = randInt(params.size());
		params[which] += randh();
		wrap(params[which], 0., 1.);
	}
	compute_scalars();

	return 0.;
}

void MyModel::compute_scalars()
{
	scalars[0] = 0.;
	scalars[1] = 0.;
	for(size_t i=0; i<params.size(); i++)
	{
		scalars[0] += -pow(params[i] - 0.5, 2);
		scalars[1] += -pow(sin(4.*M_PI*params[i]), 2);
	}
}

ostream& operator << (ostream& out, const MyModel& m)
{
	for(size_t i=0; i<m.params.size(); i++)
		out<<m.params[i]<<' ';
	return out;
}

