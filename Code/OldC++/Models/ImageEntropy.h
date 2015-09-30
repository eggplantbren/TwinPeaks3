#ifndef _ImageEntropy_
#define _ImageEntropy_

/*
* An object of this class represents a
* point in the parameter space.
*/

#include <vector>
#include <ostream>
#include "Model.h"
#include "PSF.h"

class ImageEntropy:public Model
{
	private:
		// Data and PSF
		static PSF psf;
		static std::vector< std::vector<double> > data;
		static PSF preblur;

		std::vector< std::vector<double> > image;
		void compute_scalars();	

	public:
		ImageEntropy();
		~ImageEntropy();

		void from_prior();
		double perturb();

		static void load_data();

		friend std::ostream& operator << (std::ostream& out,
							const ImageEntropy& m);
};

std::ostream& operator << (std::ostream& out, const ImageEntropy& m);

#endif

