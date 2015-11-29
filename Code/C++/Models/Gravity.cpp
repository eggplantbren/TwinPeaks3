#include "Gravity.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"

#include <cmath>

using namespace DNest3;
using namespace std;

Gravity::Gravity()
:x(100), y(100), z(100), vx(100), vy(100), vz(100)
,staleness(0)
,scalars(2)
{

}

void Gravity::from_prior(RNG& rng)
{
	for(size_t i=0; i<x.size(); i++)
	{
		x[i] = -10. + 20.*rng.rand();
		y[i] = -10. + 20.*rng.rand();
		z[i] = -10. + 20.*rng.rand();
		vx[i] = -10. + 20.*rng.rand();
		vy[i] = -10. + 20.*rng.rand();
		vz[i] = -10. + 20.*rng.rand();
	}

	refresh();
	compute_scalars();
}

void Gravity::increment(int i, int sign)
{
	KE += sign*0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	L += sign*(x[i]*vy[i] - y[i]*vx[i]);

	double rsq;
	for(size_t j=0; j<x.size(); j++)
	{
		if(i != static_cast<int>(j))
		{
			rsq = pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2)
				+ pow(z[i] - z[j], 2);
			if(rsq <= 0.01)
				rsq = 0.01;
			PE += -1./sqrt(rsq)*sign;
		}
	}
}

double Gravity::perturb(RNG& rng)
{
	int reps = 1;
	if(rng.rand() <= 0.5)
		reps += 1 + rng.rand_int(9);

	int which;
	for(int i=0; i<reps; i++)
	{
		which = rng.rand_int(x.size());

		// Remove effect of this particle
		increment(which, -1);

		x[which] += 20.*rng.randh();
		x[which] = mod(x[which] + 10., 20.) - 10.;

		y[which] += 20.*rng.randh();
		y[which] = mod(y[which] + 10., 20.) - 10.;

		z[which] += 20.*rng.randh();
		z[which] = mod(z[which] + 10., 20.) - 10.;

		vx[which] += 20.*rng.randh();
		vx[which] = mod(vx[which] + 10., 20.) - 10.;

		vy[which] += 20.*rng.randh();
		vy[which] = mod(vy[which] + 10., 20.) - 10.;

		vz[which] += 20.*rng.randh();
		vz[which] = mod(vz[which] + 10., 20.) - 10.;

		increment(which, +1);
	}

	staleness++;
	if(staleness >= 1000)
		refresh();

	compute_scalars();
	return 0.;
}

void Gravity::refresh()
{	
	KE = 0.; PE = 0.; L = 0.;
	for(size_t i=0; i<x.size(); i++)
		increment(i, 1);
	PE *= 0.5;
	staleness = 0;
}

void Gravity::compute_scalars()
{
	scalars[0] = -(KE + PE);
	scalars[1] = L;
}

void Gravity::write_text(ostream& out) const
{
	const Gravity& e = *this;
	for(size_t i=0; i<e.x.size(); i++)
		out<<e.x[i]<<' ';
	for(size_t i=0; i<e.y.size(); i++)
		out<<e.y[i]<<' ';
	for(size_t i=0; i<e.z.size(); i++)
		out<<e.z[i]<<' ';
	for(size_t i=0; i<e.vx.size(); i++)
		out<<e.vx[i]<<' ';
	for(size_t i=0; i<e.vy.size(); i++)
		out<<e.vy[i]<<' ';
	for(size_t i=0; i<e.vz.size(); i++)
		out<<e.vz[i]<<' ';

	out<<e.staleness;
}

