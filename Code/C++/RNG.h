#ifndef _RNG_
#define _RNG_

#include <random>

class RNG
{
	private:
		// 64-bit Mersenne Twister
		std::mt19937_64 twister;

		// For uniform distribution
		std::uniform_real_distribution<double> uniform;

		// For normal distribution
		std::normal_distribution<double> normal;

	public:
		// Constructor
		RNG();

		// Set the seed (obviously)
		void set_seed(unsigned int seed);

		// Uniform(0, 1)
		double rand();

		// Normal(0, 1)
		double randn();

		// My favourite heavy-tailed distribution
		double randh();

		// Integer from {0, 1, 2, ..., N-1}
		int rand_int(int N);
};

#endif

