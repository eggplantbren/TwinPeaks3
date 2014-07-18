
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles, int mcmc_steps, int thin)
:particles(num_particles)
,threshold(particles[0].get_scalars().size(), std::vector<double>(2))
,mcmc_steps(mcmc_steps)
,thin(thin)
,iteration(0)
{

}

template<class Type>
void Sampler<Type>::initialise()
{
	for(size_t i=0; i<particles.size(); i++)
	{
		particles[i].from_prior();
		particles[i].from_prior_tiebreakers();
	}

	// Set initial threshold
	for(size_t i=0; i<threshold.size(); i++)
	{
		threshold[i][0] = -1E300;
		threshold[i][1] = 0.;
	}

	// Generate direction
	direction.resize(particles[0].get_scalars().size());
	double tot = 0.;
	double sig = 5.*DNest3::randomU();
	for(size_t i=0; i<direction.size(); i++)
	{
		direction[i] = exp(sig*DNest3::randn());
		tot += direction[i];
	}
	for(size_t i=0; i<direction.size(); i++)
		direction[i] /= tot;
}

template<class Type>
void Sampler<Type>::update()
{
	// Open files
	std::fstream logw_file("logw.txt", std::ios::out|std::ios::app);
	std::fstream scalars_file("scalars.txt", std::ios::out|std::ios::app);
	std::fstream scalars_thinned_file("scalars_thinned.txt", std::ios::out|std::ios::app);
	std::fstream logw_thinned_file("logw_thinned.txt", std::ios::out|std::ios::app);
	std::fstream sample_file("sample.txt", std::ios::out|std::ios::app);

	// Choose a scalar
	int which_scalar;
	do
	{
		which_scalar = DNest3::randInt(threshold.size());
	}while(DNest3::randomU() >= direction[which_scalar]);

	// Find worst particle
	int worst = find_worst(which_scalar);

	// Write out its prior weight and scalars
	double logw = iteration*log((double)particles.size()/(particles.size()+1));
	logw_file<<logw<<std::endl;
	for(size_t i=0; i<particles[worst].get_scalars().size();  i++)
		scalars_file<<particles[worst].get_scalars()[i]<<' ';
	scalars_file<<std::endl;

	// Save to thinned files with probability 1/thin
	if(DNest3::randomU() <= 1./thin)
	{
		for(size_t i=0; i<particles[worst].get_scalars().size();  i++)
			scalars_thinned_file<<particles[worst].get_scalars()[i]<<' ';
		scalars_thinned_file<<std::endl;
		logw_thinned_file<<logw<<std::endl;
		sample_file<<particles[worst]<<std::endl;
	}

	// Close files
	logw_file.close(); scalars_file.close();
	scalars_thinned_file.close(); logw_thinned_file.close(); sample_file.close();

	// Set the new threshold
	threshold[which_scalar][0] = particles[worst].get_scalars()[which_scalar];
	threshold[which_scalar][1] = particles[worst].get_tiebreakers()[which_scalar];

	std::cout<<"# Iteration "<<(iteration+1)<<" of run with direction = (";
	for(size_t i=0; i<direction.size(); i++)
	{
		std::cout<<direction[i];
		if(i == direction.size() - 1)
			std::cout<<")."<<std::endl;
		else
			std::cout<<", ";
	}

	std::cout<<"# logw = "<<logw<<", threshold = (";
	for(size_t i=0; i<threshold.size(); i++)
	{
		std::cout<<threshold[i][0];
		if(i != (threshold.size()-1))
			std::cout<<", ";
	}
	std::cout<<")."<<std::endl;
	std::cout<<"# Evolving...";

//	// Copy a survivor
//	int which;
//	do
//	{
//		which = DNest3::randInt(particles.size());
//	}while(which == worst);
//	particles[worst] = particles[which];

	// Evolve the particle
	int accepts = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		Type proposal = particles[worst];
		double logH = proposal.perturb();
		logH += proposal.perturb_tiebreakers();
		if(proposal.is_above(threshold) &&
				DNest3::randomU() <= exp(logH))
		{
			particles[worst] = proposal;
			accepts++;
		}
	}
	std::cout<<"done. Accepted "<<accepts<<"/"<<mcmc_steps<<".";
	std::cout<<std::endl<<std::endl<<std::endl;

	iteration++;
}

template<class Type>
int Sampler<Type>::find_worst(int which_scalar) const
{
	int worst = 0;
	for(size_t i=0; i<particles.size(); i++)
	{
		bool is_worse = false;
		if(particles[i].get_scalars()[which_scalar] <
				particles[worst].get_scalars()[which_scalar])
			is_worse = true;

		if(particles[i].get_scalars()[which_scalar] ==
				particles[worst].get_scalars()[which_scalar])
			if(particles[i].get_tiebreakers()[which_scalar] <
				particles[worst].get_tiebreakers()[which_scalar])
			is_worse = true;

		if(is_worse)
			worst = i;
	}
	return worst;
}

