#include "Atoms.h"
#include "Utils.h"
#include <cmath>

using namespace std;

Atoms::Atoms()
:x(50), y(50), z(50)
,terms1(50, vector<double>(50))
,terms2(50, vector<double>(50))
,scalars(2)
{

}

Atoms::~Atoms()
{

}

void Atoms::from_prior(RNG& rng)
{
	for(size_t i=0; i<x.size(); i++)
	{
		x[i] = rng.rand();
		y[i] = rng.rand();
		z[i] = rng.rand();
	}

	for(size_t i=0; i<terms1.size(); i++)
		for(size_t j=(i+1); j<terms1[i].size(); j++)
			calculate_PE(i, j);

	calculate_PE();
	compute_scalars();
}

void Atoms::compute_scalars()
{
	scalars[0] = -PE1;
	scalars[1] = -PE2;
}

void Atoms::calculate_PE()
{
	PE1 = 0.;
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=(i+1); j<x.size(); j++)
			PE1 += terms1[i][j];

	PE2 = 0.;
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=(i+1); j<x.size(); j++)
			PE2 += terms2[i][j];
}

void Atoms::calculate_PE(int i, int j)
{
	constexpr double Rmsq = pow(0.02, 2);
	double rsq = pow(x[i] - x[j], 2) + pow(z[i] - z[j], 2);
	terms1[i][j] = pow(Rmsq/rsq, 6) - 2.*pow(Rmsq/rsq, 3);
	terms2[i][j] = rsq;
}

double Atoms::perturb(RNG& rng)
{
	int which = rng.rand_int(x.size());
	int what = rng.rand_int(3);

	if(what == 0)
	{
		x[which] += rng.randh();
		wrap(x[which], 0., 1.);
	}
	if(what == 1)
	{
		y[which] += rng.randh();
		wrap(y[which], 0., 1.);
	}
	if(what == 2)
	{
		z[which] += rng.randh();
		wrap(z[which], 0., 1.);
	}

	for(int i=0; i<which; i++)
		calculate_PE(i, which);

	for(int j=(which+1); j<int(x.size()); j++)
		calculate_PE(which, j);


	calculate_PE();
	compute_scalars();

	return 0.;
}

void Atoms::write_text(std::ostream& out) const
{
	for(size_t i=0; i<x.size(); i++)
		out<<x[i]<<' ';
	for(size_t i=0; i<y.size(); i++)
		out<<y[i]<<' ';
	for(size_t i=0; i<z.size(); i++)
		out<<z[i]<<' ';
}

