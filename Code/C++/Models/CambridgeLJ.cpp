#include "CambridgeLJ.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

CambridgeLJ::CambridgeLJ()
:scalars(2)
{

}

CambridgeLJ::~CambridgeLJ()
{

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

