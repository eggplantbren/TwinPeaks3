#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

#include <vector>
#include <list>
#include <fstream>
#include "RNG.h"
#include "ScalarType.h"

namespace TwinPeaks
{

/*
* An object of this class is a Twinpeaks sampler.
*/

template<class MyModel>
class Sampler
{
	private:
        static constexpr size_t num_scalars = 2;

        // Sampler tuning parameters
        const unsigned int num_particles;
        const unsigned int mcmc_steps;

		// The random number generators to use
		std::vector<RNG> rngs;

		/******** The particles, and metadata about them *******/
		std::vector<MyModel> particles;
		std::vector< std::vector<ScalarType> > scalars;

        // Indices that sort by the scalars (i.e. an argsort)
        std::vector< std::vector<size_t> > indices;
        std::vector< std::vector<size_t> > ranks;

        // UCCs of the particles (the non-integer part is for tiebreaking)
        std::vector<double> particle_uccs;
        /******** END particle stuff ********/

	public:
		// Constructor
		Sampler(unsigned int num_particles, unsigned int mcmc_steps,
                const std::vector<RNG>& rngs);

		// Call from_prior on all the particles
		void initialise();

		// Do an iteration of Nested Sampling
		void do_iteration();
};

} // namespace TwinPeaks

#include "SamplerImpl.h"

#endif

