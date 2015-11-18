#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <list>
#include <fstream>
#include "RNG.h"
#include "ScalarType.h"

/*
* An object of this class is a classic Nested Sampler.
*/

template<class MyModel>
class Sampler
{
	private:
		// The random number generator to use
		RNG rng;

		// The particles, and information about them
		int num_particles;
		std::vector<MyModel> particles;
		std::vector< std::vector<ScalarType> > scalars;

		// Forbidden rectangles
		std::list< std::vector<ScalarType> > rects;

		// Number of equilibration steps
		int mcmc_steps;

		// Save interval
		int save_interval;

		// Whether from_prior has been called on all the particles
		bool initialised;

		// Count number of iterations done
		int iteration;

		// Remaining prior mass
		long double log_prior_mass;

		// Do MCMC to equilibrate a particle
		int refresh_particle(int which);

		// Check s against rects
	public:
		bool is_okay(const std::vector<ScalarType>& s);

		// Remove redundant rectangles (only do this just after adding one)
		void prune_rectangles();

		// Functions to determine whether a point is within another point's
		// rectangle.
		static bool is_in_lower_rectangle(const std::vector<ScalarType>& s1,
											const std::vector<ScalarType>& s2);
		static bool is_in_upper_rectangle(const std::vector<ScalarType>& s1,
											const std::vector<ScalarType>& s2);

	public:
		// Constructor
		Sampler(const RNG& rng, int num_particles, int mcmc_steps,
														int save_interval);

		// Set RNG seed
		void set_rng_seed(unsigned int seed);

		// Call from_prior on all the particles
		void initialise();

		// Do an iteration of Nested Sampling
		void do_iteration();

		// Do an indefinite number of iterations
		void run();

		// Do a specified number of iterations
		double run(int iterations);

		// Open and close (i.e. clear the contents of)
		// the output files. If not called, it will append
		// to whatever's already in them
		void clear_output_files();

		// Setter and getter for the RNG
		void set_rng(const RNG& rng)
		{ this->rng = rng; }
		RNG get_rng() const
		{ return rng; }
};

#include "SamplerImpl.h"

#endif

