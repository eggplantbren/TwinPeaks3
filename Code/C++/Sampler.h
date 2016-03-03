#ifndef TwinPeaks_Sampler
#define TwinPeaks_Sampler

#include <vector>
#include <list>
#include <fstream>
#include "Context.h"
#include "RNG.h"
#include "ScalarType.h"

/*
* An object of this class is a TwinPeaks sampler
*/

namespace TwinPeaks
{

template<class MyModel>
class Sampler
{
	private:
		// The random number generator to use
		std::vector<RNG> rngs;

		// The particles, and information about them
		int num_particles;
		std::vector<MyModel> particles;
		std::vector<ScalarType> scalar1;
        std::vector<ScalarType> scalar2;

		std::vector<ScalarType> scalar1_sorted;
        std::vector<ScalarType> scalar2_sorted;

        // UCCS and particle uccs
        std::vector< std::vector<unsigned short> > uccs;
        std::vector<unsigned short> particle_uccs;

		// Backup
		std::vector<MyModel> backup_particles;
		std::vector< std::vector<ScalarType> > backup_scalars;

		// Forbidden rectangles
		Context context;

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



        /* Private member functions, just used to divide up the job
           of doing a TwinPeaks iteration. */

        // Calculate the uccs and particle_uccs
        void calculate_uccs();
        unsigned short choose_ucc_threshold() const;

		// Do MCMC to equilibrate a particle
		int refresh_particle(int which, int which_rng);

		// Do MCMC to equilibrate a set of particles
		void refresh_particles(const std::vector<int>& indices, int which_rng,
										int& accepts);

	public:
		bool is_okay(const std::vector<ScalarType>& s);

		// Remove redundant rectangles
		void prune_rectangles(const std::vector<ScalarType>& latest);

	public:
		// Constructor
		Sampler(const std::vector<RNG>& rngs, int num_particles, int mcmc_steps,
														int saves_per_iteration);

		// Call from_prior on all the particles
		void initialise();

		// Do an iteration of Nested Sampling
		double do_iteration();

		// Do an indefinite number of iterations
		void run();

		// Do a specified number of iterations
		double run(int iterations);

		// Open and close (i.e. clear the contents of)
		// the output files. If not called, it will append
		// to whatever's already in them
		void clear_output_files();
};

} // namespace TwinPeaks

#include "SamplerImpl.h"

#endif

