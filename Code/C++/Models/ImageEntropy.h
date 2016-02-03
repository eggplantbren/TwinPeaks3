#ifndef _ImageEntropy_
#define _ImageEntropy_

/*
* An object of this class represents a
* point in the parameter space.
*/

#include <vector>
#include <ostream>
#include "PSF.h"
#include "RNG.h"

class ImageEntropy
{
	private:
		// Data and PSF
		static PSF psf;
		static std::vector< std::vector<double> > data;
		static PSF preblur;

		std::vector< std::vector<double> > image;
		void compute_scalars();
		std::vector<double> scalars;

	public:
		ImageEntropy();

		void from_prior(RNG& rng);
		double perturb(RNG& rng);

		static void load_data();

		void write_text(std::ostream& out) const;
		const std::vector<double>& get_scalars() const
		{ return scalars; }
};

#endif

