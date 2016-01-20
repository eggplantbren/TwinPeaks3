#include "Atoms.h"
#include "Utils.h"
#include <cmath>

using namespace std;

Atoms::Atoms()
:x(num_atoms), y(num_atoms), z(num_atoms)
,terms1(num_atoms, vector<double>(num_atoms))
,terms2(num_atoms, vector<double>(num_atoms))
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
	double rsq = pow(x[i] - x[j], 2) + pow(z[i] - z[j], 2);
	terms1[i][j] = 4*(pow(1./rsq, 6) - 2.*pow(1./rsq, 3));
	terms2[i][j] = rsq;
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

