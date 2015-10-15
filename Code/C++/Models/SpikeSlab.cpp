#include "SpikeSlab.h"
#include "Utils.h"
#include <cmath>
#include <limits>

using namespace std;

SpikeSlab::SpikeSlab()
:params(20)
,scalars(2)
{

}

void SpikeSlab::from_prior(RNG& rng)
{
	for(double& x:params)
		x = rng.rand();
	compute_scalars();
}

double SpikeSlab::perturb(RNG& rng)
{
	int which, count;
	if(rng.rand() <= 0.5)
		count = 0;
	else
		count = static_cast<int>(pow(10., 2.*rng.rand()));

	for(int i=0; i<count; i++)
	{
		which = rng.rand_int(params.size());
		params[which] += rng.randh();
		wrap(params[which], 0., 1.);
	}

	compute_scalars();
	return 0.;
}

void SpikeSlab::compute_scalars()
{
	scalars[0] = 0.;
	scalars[1] = 0.;


	double u = 0.1;
	double v = 0.01;
	double C1 = -log(u) - log(2*M_PI);
	double C2 = -log(v) - log(2*M_PI);
	double C3 = log(100.);
	double uu = u*u;
	double vv = v*v;

	double temp1, temp2;
	for(const double& x: params)
	{
		temp1 = C1 - 0.5*pow(x - 0.5, 2)/uu;
		temp2 = C2 - 0.5*pow(x - 0.5, 2)/vv;
		scalars[0] += logsumexp(C3 + temp1, temp2);
	}
	scalars[1] = scalars[0];
}

void SpikeSlab::write_text(std::ostream& out) const
{
	for(size_t i=0; i<params.size(); i++)
		out<<params[i]<<' ';
}

