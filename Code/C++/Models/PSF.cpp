#include "PSF.h"
#include "Utils.h"
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;
using namespace arma;

PSF::PSF(int size)
:size(size)
,pixels(size, vector<double>(size, 0.))
,fft_ready(false)
{
	assert(size%2 == 1);
	pixels[size/2][size/2] = 1.;
}

void PSF::set_size(int new_size)
{
	size = new_size;
	pixels.assign(size, vector<double>(size, 0.));
	pixels[size/2][size/2] = 1.;
}

void PSF::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# ERROR: couldn't open file "<<filename<<"."<<endl;
		return;
	}
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			fin>>pixels[i][j];
	fin.close();
	normalise();
}

void PSF::normalise()
{
	double sum = 0.;
	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			sum += pixels[i][j];

	for(int i=0; i<size; i++)
		for(int j=0; j<size; j++)
			pixels[i][j] /= sum;
}

void PSF::calculate_fft(int Ni, int Nj)
{
	// Make the psf the same size as the image
	mat psf(Ni, Nj);
	psf.zeros();

	int ni = pixels.size();
	int nj = pixels[0].size();

	int m, n;
	for(int i=0; i<ni; i++)
	{
		m = mod(i - ni/2, Ni);
		for(int j=0; j<nj; j++)
		{
			n = mod(j - nj/2, Nj);
			psf(m, n) = pixels[i][j];
		}
	}

	fft_of_psf = fft2(psf);
	fft_ready = true;
}

void PSF::blur_image(vector< vector<double> >& img) const
{
	// Make result image. Assume img is rectangular...
	vector< vector<double> > result(img.size(),
					vector<double>(img[0].size(), 0.));

	int h = size/2;
	int ii, jj;
	int M = static_cast<int>(img.size());
	int N = static_cast<int>(img[0].size());

	for(int i=0; i<M; i++)
	{
		for(int j=0; j<N; j++)
		{
			if(img[i][j] != 0.)
			{
				for(int m=0; m<size; m++)
				{
					ii = i + m - h;
					for(int n=0; n<size; n++)
					{
						jj = j + n - h;
						if(ii >= 0 && ii < M &&
							jj >= 0 && jj < N)
							result[ii][jj] +=
							img[i][j]*pixels[m][n];
					}
				}
			}
		}
	}

	img = result;
}

void PSF::blur_image2(vector< vector<double> >& img) const
{
	if(!fft_ready)
		cerr<<"# Blurring failed."<<endl;

	// Copy the image into an Armadillo matrix
	mat A(img.size(), img[0].size());
	for(size_t i=0; i<img.size(); i++)	
		for(size_t j=0; j<img[0].size(); j++)
			A(i, j) = img[i][j];

	// Do the fft of it
	cx_mat B = fft2(A);

	// Multiply the two ffts
	for(size_t i=0; i<img.size(); i++)
		for(size_t j=0; j<img[0].size(); j++)
			B(i, j) *= fft_of_psf(i, j);

	// Do the inverse fft
	B = ifft2(B);

	// Put back in img
	for(size_t i=0; i<img.size(); i++)
		for(size_t j=0; j<img[0].size(); j++)
			img[i][j] = real(B(i, j));

}


void PSF::test()
{
	PSF psf(5);
	psf.load("psf.txt");
	psf.calculate_fft(20, 20);

	// Point image
	vector< vector<double> > pixels(20, vector<double>(20, 0.));
	pixels[10][10] = 1.;

	psf.blur_image(pixels);
	// Print image and reset it to zero
	for(size_t i=0; i<pixels.size(); i++)
	{
		for(size_t j=0; j<pixels.size(); j++)
		{
			cout<<pixels[i][j]<<' ';
			pixels[i][j] = 0.;
		}
	}
	cout<<endl;

	// Do it again with ffts
	pixels[10][10] = 1.;

	psf.blur_image2(pixels);
	for(size_t i=0; i<pixels.size(); i++)
		for(size_t j=0; j<pixels.size(); j++)
			cout<<pixels[i][j]<<' ';
	cout<<endl;
}


