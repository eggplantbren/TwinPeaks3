#include "CambridgeLJ.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>
#include <cassert>

using namespace std;
using namespace DNest3;

// To call Rob's Fortran functions
extern "C"
{
	void __twin_peaks_routines_MOD_initialise_config(double* s, double* h0,
					double* V, double* ener, double* min_height,
					double* Vmax, int* N, double* cutoff,
					double* univar_3Nplus10001, double* univar_60k,
					double* univar_30k, double* nvar_60k, bool* flat_v_prior);
}

const int CambridgeLJ::N = 100;

CambridgeLJ::CambridgeLJ()
:scalars(2)
{
	s = new double[3*N];
}

CambridgeLJ::CambridgeLJ(const CambridgeLJ& other)
{
	s = new double[3*N];
	*this = other;
}

CambridgeLJ& CambridgeLJ::operator = (const CambridgeLJ& other)
{
	for(int i=0; i<N; i++)
		s[i] = other.s[i];
	scalars = other.scalars;
	return *this;
}

CambridgeLJ::~CambridgeLJ()
{
	delete[] s;
}

void CambridgeLJ::from_prior(RNG& rng)
{
	compute_scalars();
}

void CambridgeLJ::compute_scalars()
{
	scalars[0] = 0.;
	scalars[1] = 0.;
}

double CambridgeLJ::perturb(RNG& rng)
{
	double logH = 0.;

	compute_scalars();
	return 0.;
}

void CambridgeLJ::write_text(std::ostream& out) const
{

}

