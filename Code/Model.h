#ifndef _Model_
#define _Model_

#include <vector>

/*
* Abstract class defining what is needed to implement a model.
* This will be like DNest. An object of a subclass represents a
* Model in parameter space.
*/

class Model
{
	protected:
		int num_scalars;
		std::vector<double> scalars;
		std::vector<double> tiebreakers;

	public:
		// Pass in the number of scalars
		Model(int num_scalars);

		virtual ~Model();

		// Generate the Model from the prior
		virtual void from_prior() = 0;

		// Do a Metropolis-Hastings proposal that satisfies detailed
		// balance wrt the prior. Return log of any ratios needed to
		// achieve this.
		virtual double perturb() = 0;

		// Methods for tiebreakers
		void from_prior_tiebreakers();
		double perturb_tiebreakers();

		// Are all scalars greater than or equal to the threshold
		// values?
		bool is_above(const std::vector< std::vector<double> >&
					threshold) const;

		// Getters for scalars and tiebreakers
		int get_num_scalars() const { return num_scalars; }
		const std::vector<double>& get_scalars() const
		{
			return scalars;
		}
		const std::vector<double>& get_tiebreakers() const
		{
			return tiebreakers;
		}
};

#endif

