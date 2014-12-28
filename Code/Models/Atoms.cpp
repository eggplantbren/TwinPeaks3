#include "Atoms.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

Atoms::Atoms()
:Model(2), x(200), y(200)
{

}

Atoms::~Atoms()
{

}

void Atoms::from_prior()
{
	for(size_t i=0; i<x.size(); i++)
	{
		x[i] = randomU();
		y[i] = randomU();
	}

	calculate_PE();
	compute_scalars();
}

void Atoms::compute_scalars()
{
	scalars[0] = -PE;
	scalars[1] = 0.;
	for(size_t i=0; i<x.size(); i++)
		scalars[1] += -y[i];
//	scalars[1] /= x.size();
}

void Atoms::calculate_PE()
{
	PE = 0.;
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=(i+1); j<x.size(); j++)
			PE += calculate_PE(i, j);
}

double Atoms::calculate_PE(int i, int j)
{
	double Rmsq = pow(0.01, 2);
	double rsq = pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2);
	return pow(Rmsq/rsq, 6) - 2.*pow(Rmsq/rsq, 3);
}

double Atoms::perturb()
{
	int which = randInt(x.size());

	double diff = 0.;
	for(size_t i=0; i<x.size(); i++)
		if((int)i != which)
			diff -= calculate_PE(i, which);

	x[which] += randh();
	y[which] += randh();
	wrap(x[which], 0., 1.);
	wrap(y[which], 0., 1.);

	for(size_t i=0; i<x.size(); i++)
		if((int)i != which)
			diff += calculate_PE(i, which);

	// If fractional change is big
	if(abs(diff)/abs(PE) > 1. || randomU() < 0.01)
		calculate_PE();
	else
		PE += diff;

	compute_scalars();

	return 0.;
}

std::ostream& operator << (std::ostream& out, const Atoms& a)
{
	for(size_t i=0; i<a.x.size(); i++)
		out<<a.x[i]<<' ';
	for(size_t i=0; i<a.y.size(); i++)
		out<<a.y[i]<<' ';
	return out;
}

