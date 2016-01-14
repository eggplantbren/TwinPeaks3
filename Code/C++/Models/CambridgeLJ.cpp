#include "CambridgeLJ.h"
#include "Utils.h"
#include <cmath>
#include <cassert>

using namespace std;

// To call Rob's Fortran functions
extern "C"
{
	void __twin_peaks_routines_MOD_initialise_config(double* s, double* h0,
					double* V, double* ener, double* min_height,
					double* Vmax, int* N, double* cutoff,
					double* univar_3Nplus10001, double* univar_60k,
					double* univar_30k, double* nvar_60k, bool* flat_v_prior);
}

CambridgeLJ::CambridgeLJ()
:N(100)
,s(0., 3*N)
,scalars(2)
{
}

void CambridgeLJ::from_prior(RNG& rng)
{
	// "Outputs" (along with s)
	std::valarray<double> h0(9);
	double V, ener;

	// "constants"
	double min_height = 0.1;
	double Vmax = 10.;
	double cutoff = 20.;
	bool flat_v_prior = false;

	// Random numbers
	std::vector<double> univar_3Nplus10001(3*N+10001);
	std::vector<double> univar_60k(60000);
	std::vector<double> univar_30k(30000);
	std::vector<double> nvar_60k(60000);
	for(double& u: univar_3Nplus10001)
		u = rng.rand();
	for(double& u: univar_60k)
		u = rng.rand();
	for(double& u: univar_30k)
		u = rng.rand();
	for(double& n: nvar_60k)
		n = rng.rand();

	// Call Rob's code
	__twin_peaks_routines_MOD_initialise_config(&s[0], &h0[0], &V,
					&ener, &min_height, &Vmax, &N, &cutoff,
					&univar_3Nplus10001[0], &univar_60k[0],
					&univar_30k[0], &nvar_60k[0], &flat_v_prior);

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
	logH += rng.rand();

	compute_scalars();
	return 0.;
}

void CambridgeLJ::write_text(std::ostream& out) const
{
	out<<' ';
}

