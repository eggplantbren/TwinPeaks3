#include "ImageEntropy.h"
#include "Utils.h"
#include <cmath>
#include <limits>

using namespace std;

vector< vector<double> > ImageEntropy::data;
PSF ImageEntropy::psf(3);

ImageEntropy::ImageEntropy()
:scalars(2)
{

}

void ImageEntropy::from_prior(RNG& rng)
{
	compute_scalars();
}

double ImageEntropy::perturb(RNG& rng)
{
	compute_scalars();
	return 0.;
}

void ImageEntropy::compute_scalars()
{
	scalars[0] = 0.;
	scalars[1] = 0.;
}

void ImageEntropy::write_text(std::ostream& out) const
{
}

