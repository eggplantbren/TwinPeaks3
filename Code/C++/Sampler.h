#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

#include <vector>
#include <list>
#include <fstream>
#include "Rectangle.h"
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
        // Sampler tuning parameters
        const unsigned int num_particles;
        const unsigned int mcmc_steps;

        // Iteration counter
        unsigned int iteration;

		// The random number generators to use
		std::vector<RNG> rngs;

		/**** The particles, and metadata about them ****/
		std::vector<MyModel> particles;
		std::vector< std::vector<ScalarType> > scalars;

        // Indices that sort by the scalars (i.e. an argsort)
        std::vector< std::vector<size_t> > indices;
        std::vector< std::vector<size_t> > ranks;

        // UCCs of the particles (the non-integer part is for tiebreaking)
        std::vector<unsigned short> particle_uccs;
        std::vector<double> particle_ucc_tiebreakers;

        /**** The forbidden rectangles ****/
        std::list<Rectangle> forbidden_rectangles;

        /**** Private member functions ****/
        void calculate_uccs();
        void forbid_rectangle(size_t which_particle, bool unique);
        void replace_particle(size_t which_particle);
        double log_prob(const std::vector<ScalarType>& s);

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

