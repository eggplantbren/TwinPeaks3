#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <fstream>
#include "RNG.h"

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

		// Forbidden rectangles
		std::vector< std::vector<double> > rects;

		// Number of equilibration steps
		int mcmc_steps;

		// Whether from_prior has been called on all the particles
		bool initialised;

		// Count number of iterations done
		int iteration;

		// Do MCMC to equilibrate a particle
		void refresh_particle(int which);

		// Check s against rects
		bool is_okay(const std::vector<double>& s);

		// Function to determine whether a point is within another point's
		// rectangle.
		static bool is_in_rectangle(const std::vector<double>& s1,
									const std::vector<double>& s2);

	public:
		// Constructor
		Sampler(const RNG& rng, int num_particles, int mcmc_steps);

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

