#ifndef _PSF_
#define _PSF_

#include <vector>
#include <armadillo>

class PSF
{
	private:
		int size;
		std::vector< std::vector<double> > pixels;

		// FFT of the PSF
		arma::cx_mat fft_of_psf;
		bool fft_ready;

	public:
		PSF(int size);
		// Uses loopy method
		void blur_image(std::vector< std::vector<double> >& img) const;
		// Uses fft method
		void blur_image2
			(std::vector< std::vector<double> >& img) const;
		void normalise();
		void calculate_fft(int Ni, int Nj);
		void load(const char* filename);
		void set_size(int new_size);

		// Getter for pixels
		const std::vector< std::vector<double> >& get_pixels() const
		{ return pixels; }

		// Unit test
		static void test();
};

#endif

