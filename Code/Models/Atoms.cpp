#include "Atoms.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

Atoms::Atoms()
:Model(2), x(200), y(200), z(200)
,terms(200, vector<double>(200))
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
		z[i] = randomU();
	}

	for(size_t i=0; i<terms.size(); i++)
		for(size_t j=(i+1); j<terms[i].size(); j++)
			calculate_PE(i, j);

	calculate_PE();
	compute_scalars();
}

void Atoms::compute_scalars()
{
	scalars[0] = -PE;

	double Lx, Ly, Lz;
	Lx = *max_element(x.begin(), x.end());
	Ly = *max_element(y.begin(), y.end());
	Lz = *max_element(z.begin(), z.end());

	scalars[1] = -log(Lx*Ly*Lz);
}

void Atoms::calculate_PE()
{
	PE = 0.;
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=(i+1); j<x.size(); j++)
			PE += terms[i][j];
}

void Atoms::calculate_PE(int i, int j)
{
	double Rmsq = pow(0.05, 2);
	double rsq = pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2) + pow(z[i] - z[j], 2);
	terms[i][j] = pow(Rmsq/rsq, 6) - 2.*pow(Rmsq/rsq, 3);
}

double Atoms::perturb()
{
	int which = randInt(x.size());

	x[which] += randh();
	y[which] += randh();
	z[which] += randh();
	wrap(x[which], 0., 1.);
	wrap(y[which], 0., 1.);
	wrap(z[which], 0., 1.);

	for(int i=0; i<which; i++)
		calculate_PE(i, which);

	for(int j=(which+1); j<int(x.size()); j++)
		calculate_PE(which, j);


	calculate_PE();
	compute_scalars();

	return 0.;
}

std::ostream& operator << (std::ostream& out, const Atoms& a)
{
	for(size_t i=0; i<a.x.size(); i++)
		out<<a.x[i]<<' ';
	for(size_t i=0; i<a.y.size(); i++)
		out<<a.y[i]<<' ';
	for(size_t i=0; i<a.z.size(); i++)
		out<<a.z[i]<<' ';
	return out;
}

