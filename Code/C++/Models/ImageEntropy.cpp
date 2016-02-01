#include "ImageEntropy.h"
#include "Utils.h"
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

vector< vector<double> > ImageEntropy::data(100, vector<double>(100));
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

void ImageEntropy::load_data(const char* filename)
{
	fstream fin(filename, ios::in);
	for(size_t i=0; i<data.size(); ++i)
		for(size_t j=0; j<data[i].size(); ++j)
			fin>>data[i][j];
	fin.close();
}

void ImageEntropy::load_psf(const char* filename)
{
	psf.load(filename);
}

