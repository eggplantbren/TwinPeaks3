#ifndef _Sampler_
#define _Sampler_

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

template<class Type>
class Sampler
{
	private:
		std::vector<Type> particles;
		std::vector<double> threshold;
		int mcmc_steps;
		int iteration;

		std::fstream logw_file, scalars_file;

		int find_worst() const;

	public:
		Sampler(int num_particles, int mcmc_steps);

		void initialise();
		void update();

};


/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles, int mcmc_steps)
:particles(num_particles)
,mcmc_steps(mcmc_steps)
,iteration(0)
,logw_file()
,scalars_file()
{

}

template<class Type>
void Sampler<Type>::initialise()
{
	for(size_t i=0; i<particles.size(); i++)
		particles[i].from_prior();

	threshold.assign(particles[0].get_scalars().size(), -1E300);
}

template<class Type>
void Sampler<Type>::update()
{
	// Open files
	logw_file.open("logw.txt",
		(iteration == 0)?(std::ios::out):(std::ios::out|std::ios::app));
	scalars_file.open("scalars.txt",
		(iteration == 0)?(std::ios::out):(std::ios::out|std::ios::app));

	// Find worst particle
	int worst = find_worst();

	// Write out its prior weight and scalars
	double logw = iteration*log((double)particles.size()/(particles.size()+1));
	logw_file<<logw<<std::endl;
	for(size_t i=0; i<particles[worst].get_scalars().size();  i++)
		scalars_file<<particles[worst].get_scalars()[i]<<' ';
	scalars_file<<std::endl;

	// Close files
	logw_file.close(); scalars_file.close();

	// Set the new threshold
	threshold[0] = particles[worst].get_scalars()[0];

	std::cout<<"# Iteration "<<(iteration+1)<<". Threshold = ";
	std::cout<<threshold[0]<<". Evolving...";

	// Evolve the particle
	int accepts = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		Type proposal = particles[worst];
		double logH = proposal.perturb();
		if(proposal.is_above(threshold) &&
				DNest3::randomU() <= exp(logH))
		{
			particles[worst] = proposal;
			accepts++;
		}
	}
	std::cout<<"done. Accepted "<<accepts<<"/"<<mcmc_steps<<"."<<std::endl;

	iteration++;
}

template<class Type>
int Sampler<Type>::find_worst() const
{
	int worst = 0;
	for(size_t i=0; i<particles.size(); i++)
	{
		if(particles[i].get_scalars()[0] <
				particles[worst].get_scalars()[0])
			worst = i;
	}
	return worst;
}

#endif

