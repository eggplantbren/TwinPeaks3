#include "Atoms.h"
#include "Utils.h"
#include <cmath>

using namespace std;

Atoms::Atoms()
:x(num_atoms), y(num_atoms), z(num_atoms)
,terms1(num_atoms, vector<double>(num_atoms))
,scalars(2)
{

}

void Atoms::from_prior(RNG& rng)
{
	for(size_t i=0; i<x.size(); i++)
	{
		x[i] = L*rng.rand();
		y[i] = L*rng.rand();
		z[i] = L*rng.rand();
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
	scalars[1] = 0.;

	// Polymer potential
	double k = 36*pow(2., 2./3);
	double c = pow(2., 1./6);
	double r;
	// First sum
	for(int i=0; i<(num_atoms-1); ++i)
	{
		r = sqrt(pow(x[i] - x[i+1], 2) + pow(y[i] - y[i+1], 2)
						+ pow(z[i] - z[i+1], 2));
		scalars[1] -= 0.5*k*pow(r - c, 2);
	}

	// Second sum
	double dotprod;
	for(int i=1; i<(num_atoms-1); ++i)
	{
		dotprod = (x[i] - x[i-1])*(x[i+1] - x[i])
					+ (y[i] - y[i-1])*(y[i+1] - y[i])
					+ (z[i] - z[i-1])*(z[i+1] - z[i]);
		scalars[1] -= 0.5*k*pow(dotprod - 1., 2);
	}
}

void Atoms::calculate_PE()
{
	PE1 = 0.;
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=(i+1); j<x.size(); j++)
			PE1 += terms1[i][j];
}

void Atoms::calculate_PE(int i, int j)
{
	double rsq = pow(x[i] - x[j], 2) + pow(z[i] - z[j], 2);
	terms1[i][j] = 4*(pow(1./rsq, 6) - 2.*pow(1./rsq, 3));
}

double Atoms::perturb(RNG& rng)
{
	int which = rng.rand_int(x.size());
	int what = rng.rand_int(3);

	if(what == 0)
	{
		x[which] += L*rng.randh();
		wrap(x[which], 0., L);
	}
	if(what == 1)
	{
		y[which] += L*rng.randh();
		wrap(y[which], 0., L);
	}
	if(what == 2)
	{
		z[which] += L*rng.randh();
		wrap(z[which], 0., L);
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

