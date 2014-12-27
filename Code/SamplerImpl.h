
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles)
:num_particles(num_particles)
,particles(num_particles)
,log_prior_mass(0.)
{

}

template<class Type>
int Sampler<Type>::badness(const Type& particle) const
{
	int count = 0;
	for(size_t i=0; i<thresholds.size(); i++)
		if(is_below(particle.get_scalars(), thresholds[i]))
			count++;
	return count;
}

template<class Type>
void Sampler<Type>::initialise()
{
	for(size_t i=0; i<particles.size(); i++)
	{
		particles[i].from_prior();
		particles[i].from_prior_tiebreakers();
	}

	thresholds.clear();

	// Open output file
	std::fstream fout("output.txt", std::ios::out);
	fout.close();
}

// Exploration just to ensure everything's above all thresholds
template<class Type>
void Sampler<Type>::refresh()
{
	const int steps = 10000;

	std::vector<int> bad(num_particles);
	for(int i=0; i<num_particles; i++)
		bad[i] = badness(particles[i]);

	for(int i=0; i<steps; i++)
	{
		int which = DNest3::randInt(num_particles);
		Type proposal = particles[which];
		double logH = proposal.perturb();
		int proposal_badness = badness(proposal);

		if(proposal_badness <= bad[which] &&
				DNest3::randomU() <= exp(logH))
		{
			particles[which] = proposal;
			bad[which] = proposal_badness;
		}
	}

	// Make sure they're not all bad
	bool all_bad = true;
	for(int i=0; i<num_particles; i++)
		if(bad[i] == 0)
			all_bad = false;

	if(all_bad)
		std::cerr<<"# WARNING: All particles are bad."<<std::endl;

	// Resample any bad points by copying good ones
	int copy;
	for(int i=0; i<num_particles; i++)
	{
		if(bad[i] > 0)
		{
			do
			{
				copy = DNest3::randInt(num_particles);
			}while(bad[copy] > 0);
			particles[i] = particles[copy];
		}
	}
}


template<class Type>
void Sampler<Type>::explore()
{
	const int steps = 50000;
	const int skip = 10;
	std::vector< std::vector<double> > keep(steps/skip);

	for(int i=0; i<steps; i++)
	{
		int which = DNest3::randInt(num_particles);

		Type proposal = particles[which];
		double logH = proposal.perturb();
		int proposal_badness = badness(proposal);

		if(DNest3::randomU() <= exp(logH) && proposal_badness == 0)
			particles[which] = proposal;

		if(i%skip == 0)
			keep[i/skip] = particles[which].get_scalars();
	}

	create_threshold(keep);
}

template<class Type>
bool Sampler<Type>::is_below(const std::vector<double>& s1,
				const std::vector<double>& s2) const
{
	for(size_t i=0; i<s1.size(); i++)
		if(s1[i] >= s2[i])
			return false;
	return true;
}

template<class Type>
void Sampler<Type>::create_threshold(const std::vector< std::vector<double> >&
						keep)
{
	int which = 0;		// Closest element to 1-exp(-1)
	double diff = 1E300;	// Distance from 0.05

	std::vector<double> frac_below(keep.size());
	for(size_t i=0; i<keep.size(); i++)
	{
		frac_below[i] = 0.;
		for(size_t j=0; j<keep.size(); j++)
		{
			if(i != j)
				frac_below[i] += is_below(keep[j], keep[i]);
		}
		frac_below[i] /= (keep.size() - 1);
		if(fabs(frac_below[i] - 0.05) < diff)
		{
			which = i;
			diff = fabs(frac_below[i] - 0.05);
		}
	}

	std::cout<<"# New threshold = ";
	for(size_t i=0; i<keep[which].size(); i++)
		std::cout<<keep[which][i]<<' ';
	std::cout<<"."<<std::endl;
	thresholds.push_back(keep[which]);

	// Write out dead points
	double log_dead_mass = log(frac_below[which]) + log_prior_mass;
	std::fstream fout("output.txt", std::ios::out|std::ios::app);
	for(size_t i=0; i<keep.size(); i++)
	{
		if(int(i) != which && is_below(keep[i], keep[which]))
		{
			fout<<(log_dead_mass - log(keep.size()*frac_below[which]))<<' ';
			for(size_t j=0; j<keep[i].size(); j++)
				fout<<keep[i][j]<<' ';
			fout<<std::endl;
		}
	}
	fout.close();

	log_prior_mass = DNest3::logdiffexp(log_prior_mass, log_dead_mass);
	std::cout<<"# Peeling away "<<frac_below[which]<<" of the remaining prior mass."<<std::endl;
	std::cout<<"# log(remaining prior mass) = "<<log_prior_mass<<std::endl;
	std::cout<<std::endl;
}

