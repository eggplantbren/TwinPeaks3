#include "SimpleExample.h"
#include "Utils.h"
#include <cmath>
#include <limits>

using namespace std;

SimpleExample::SimpleExample()
:params(100)
,scalars(2)
{

}

void SimpleExample::from_prior(RNG& rng)
{
	for(double& x:params)
		x = rng.rand();
	compute_scalars();
}

double SimpleExample::perturb(RNG& rng)
{
	int which = rng.rand_int(params.size());
	params[which] += rng.randh();
	wrap(params[which], 0., 1.);
	compute_scalars();
	return 0.;
}

void SimpleExample::compute_scalars()
{
	scalars[0] = 0.;
	scalars[1] = 0.;
	for(const double& x: params)
	{
		scalars[0] += -x*x;
		scalars[1] += pow(sin(4.*M_PI*x), 2);
	}
}

void SimpleExample::write_text(std::ostream& out) const
{
	for(size_t i=0; i<params.size(); i++)
		out<<params[i]<<' ';
}

