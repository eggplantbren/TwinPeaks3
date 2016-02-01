#include "ImageEntropy.h"
#include "Utils.h"
#include <cmath>
#include <limits>
#include <fstream>

using namespace std;

vector< vector<double> > ImageEntropy::data(100, vector<double>(100));
PSF ImageEntropy::psf(21);

ImageEntropy::ImageEntropy()
:scalars(2)
{
}

void ImageEntropy::from_prior(RNG& rng)
{
	for(size_t i=0; i<image.size(); ++i)
		for(size_t j=0; j<image.size(); ++j)
			image[i][j] = 1E3*rng.rand();
	compute_scalars();
}

double ImageEntropy::perturb(RNG& rng)
{
	int i, j;
	int reps = static_cast<int>(pow(10., 2.*rng.rand()));
	for(int rep=0; rep<reps; ++rep)
	{
		i = rng.rand_int(image.size());
		j = rng.rand_int(image[0].size());
		image[i][j] += 1E3*rng.randh();
		wrap(image[i][j], 0., 1E3);
	}
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
	for(size_t i=0; i<image.size(); ++i)
		for(size_t j=0; j<image.size(); ++j)
			out<<image[i][j]<<' ';
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

