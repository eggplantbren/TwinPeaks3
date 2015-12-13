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
		std::vector<RNG> rngs;

		// The particles, and information about them
		int num_particles;
		std::vector<MyModel> particles;
		std::vector< std::vector<ScalarType> > scalars;

		// Backup
		std::vector<MyModel> backup_particles;
		std::vector< std::vector<ScalarType> > backup_scalars;

		// Dying flags
		std::vector<bool> dying;

		// Forbidden rectangles
		std::list< std::vector<ScalarType> > rects;

		// Number of equilibration steps
		int mcmc_steps;

		// How many samples (to save to sample.txt) per iteration
		int saves_per_iteration;

		// Whether from_prior has been called on all the particles
		bool initialised;

		// Count number of iterations done
		int iteration;

		// Remaining prior mass
		long double log_prior_mass;

		// Do MCMC to equilibrate a particle
		int refresh_particle(int which, int which_rng);

		// Do MCMC to equilibrate a set of particles
		void refresh_particles(const std::vector<int>& indices, int which_rng,
										int& accepts);

	public:
		bool is_okay(const std::vector<ScalarType>& s);

		// Remove redundant rectangles
		void prune_rectangles(const std::vector<ScalarType>& latest);

		// Functions to determine whether a point is within another point's
		// rectangle. Strict version
		static bool is_in_lower_rectangle(const std::vector<ScalarType>& s1,
											const std::vector<ScalarType>& s2);

		// Non-strict version
		static bool is_in_lower_rectangle2(const std::vector<ScalarType>& s1,
											const std::vector<ScalarType>& s2);



	public:
		// Constructor
		Sampler(const std::vector<RNG>& rngs, int num_particles, int mcmc_steps,
														int saves_per_iteration);

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
};

#include "SamplerImpl.h"

#endif

