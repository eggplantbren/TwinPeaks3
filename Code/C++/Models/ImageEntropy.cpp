#include "ImageEntropy.h"
#include <cassert>
#include <cmath>
#include "Utils.h"
#include <fstream>

using namespace std;

vector< vector<double> > ImageEntropy::data(100, vector<double>(100));
PSF ImageEntropy::psf(21);
PSF ImageEntropy::preblur(5);

void ImageEntropy::load_data()
{
	fstream fin("data.txt", ios::in);
	for(size_t i=0; i<data.size(); i++)
		for(size_t j=0; j<data[i].size(); j++)
			fin>>data[i][j];
	fin.close();

	psf.load("Models/psf.txt");
	psf.calculate_fft(100, 100);

	preblur.load("Models/preblur.txt");
	preblur.calculate_fft(100, 100);
}

ImageEntropy::ImageEntropy()
:image(100, vector<double>(100))
,scalars(2)
{

}

void ImageEntropy::from_prior(RNG& rng)
{
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			image[i][j] = rng.rand();
	compute_scalars();
}

double ImageEntropy::perturb(RNG& rng)
{
	int reps = 1;
	if(rng.rand() <= 0.9)
		reps = (int)pow(10., 3.*rng.rand());

	int ii, jj;
	for(int i=0; i<reps; ++i)
	{
		ii = rng.rand_int(image.size());
		jj = rng.rand_int(image[ii].size());
		image[ii][jj] += rng.randh();
		wrap(image[ii][jj], 0., 1.);
	}

	compute_scalars();
	return 0.;
}

void ImageEntropy::compute_scalars()
{
	// Find image total (use entropy of normalised image)
	double tot = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			tot += image[i][j];

	// Find image entropy
	double S = 0.;
	for(size_t i=0; i<image.size(); ++i)
		for(size_t j=0; j<image[i].size(); ++j)
			S += -(image[i][j]/tot)*log(image[i][j]/tot + 1E-300);
	scalars[0] = S;

	vector< vector<double> > blurred = image;
	preblur.blur_image2(blurred);
	psf.blur_image2(blurred);

	// Log likelihood
	double logL = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			logL += -0.5*pow((data[i][j] - blurred[i][j])/0.2, 2);
	scalars[0] = S;
	scalars[1] = logL;
}

void ImageEntropy::write_text(ostream& out) const
{
	vector< vector<double> > blurred = image;
	preblur.blur_image2(blurred);

	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			out<<blurred[i][j]<<' ';

	psf.blur_image2(blurred);

	// Log likelihood
	double logL = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			logL += -0.5*pow((ImageEntropy::data[i][j] - blurred[i][j])/0.2, 2);

	out<<logL<<' ';
}

