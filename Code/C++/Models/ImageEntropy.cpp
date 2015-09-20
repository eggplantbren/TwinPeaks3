#include "ImageEntropy.h"
#include <cassert>
#include <cmath>
#include <RandomNumberGenerator.h>
#include <Utils.h>
#include <fstream>

using namespace std;
using namespace DNest3;

vector< vector<double> > ImageEntropy::data(100, vector<double>(100));
PSF ImageEntropy::psf(21);

void ImageEntropy::load_data()
{
	fstream fin("data.txt", ios::in);
	for(size_t i=0; i<data.size(); i++)
		for(size_t j=0; j<data[i].size(); j++)
			fin>>data[i][j];
	fin.close();

	psf.load("psf.txt");
	psf.calculate_fft(100, 100);
}

ImageEntropy::ImageEntropy()
:Model(2)
,image(100, vector<double>(100))
{

}

ImageEntropy::~ImageEntropy()
{

}

void ImageEntropy::from_prior()
{
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			image[i][j] = randomU();
	compute_scalars();
}

double ImageEntropy::perturb()
{
	int reps = 1;
	if(randomU() <= 0.9)
		reps = (int)pow(10., 4.*randomU());

	int ii, jj;
	for(int i=0; i<reps; i++)
	{
		ii = randInt(image.size());
		jj = randInt(image[ii].size());
		image[ii][jj] += randh();
		wrap(image[ii][jj], 0., 1.);
	}

	compute_scalars();
	return 0.;
}

void ImageEntropy::compute_scalars()
{
	vector< vector<double> > blurred = image;
	psf.blur_image2(blurred);

	// Find image total (use entropy of normalised image)
	double tot = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			tot += image[i][j];

	// Find image entropy
	double S = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			S += -(image[i][j]/tot)*log(image[i][j]/tot + 1E-300);
	scalars[0] = S;

	// Log likelihood
	double logL = 0.;
	for(size_t i=0; i<image.size(); i++)
		for(size_t j=0; j<image[i].size(); j++)
			logL += -0.5*pow((data[i][j] - blurred[i][j])/0.2, 2);
	scalars[1] = 1000*S + logL;
	scalars[0] = 1000*S + logL;
}

ostream& operator << (ostream& out, const ImageEntropy& m)
{
	for(size_t i=0; i<m.image.size(); i++)
		for(size_t j=0; j<m.image[i].size(); j++)
			out<<m.image[i][j]<<' ';
	return out;
}

