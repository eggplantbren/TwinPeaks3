
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles, int steps, double peel_factor)
:num_particles(num_particles)
,steps(steps)
,peel_factor(peel_factor)
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
	std::cout<<"# Evolving particles..."<<std::flush;
	int accepts = 0;

	std::vector<int> bad(num_particles);
	std::vector<int> need_refresh;
	for(int i=0; i<num_particles; i++)
	{
		bad[i] = badness(particles[i]);
		if(bad[i] != 0)
			need_refresh.push_back(i);
	}
	if(need_refresh.size() == 0)
		return;

	int which, proposal_badness;
	Type proposal;
	double logH;
	for(int i=0; i<steps; i++)
	{
		which = need_refresh[DNest3::randInt(need_refresh.size())];

		proposal = particles[which];
		logH = proposal.perturb();
		proposal_badness = badness(proposal);

		if(proposal_badness <= bad[which] &&
				DNest3::randomU() <= exp(logH))
		{
			particles[which] = proposal;
			bad[which] = proposal_badness;
			accepts++;
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
	need_refresh.clear();
	for(int i=0; i<num_particles; i++)
	{
		if(bad[i] > 0)
		{
			do
			{
				copy = DNest3::randInt(num_particles);
			}while(bad[copy] > 0);
			particles[i] = particles[copy];
			bad[i] = bad[copy];
			need_refresh.push_back(i);
		}
	}

	// Evolve any copied particles
	for(int i=0; i<steps; i++)
	{
		which = need_refresh[DNest3::randInt(need_refresh.size())];

		proposal = particles[which];
		logH = proposal.perturb();
		proposal_badness = badness(proposal);

		if(proposal_badness <= bad[which] &&
				DNest3::randomU() <= exp(logH))
		{
			particles[which] = proposal;
			bad[which] = proposal_badness;
		}
		accepts++;
	}
	std::cout<<"done. Acceptance rate = "<<accepts<<"/"<<(2*steps)<<"."<<std::endl<<std::endl;
}

template<class Type>
void Sampler<Type>::remove_redundant_thresholds()
{
	top:
	for(size_t i=0; i<thresholds.size()-1; i++)
	{
		if(is_below(thresholds[i], thresholds.back()))
		{
			thresholds.erase(thresholds.begin() + i);
			goto top;
		}
	}
}

template<class Type>
void Sampler<Type>::explore()
{
	std::vector< std::vector<double> > keep(num_particles);

	for(int i=0; i<num_particles; i++)
		keep[i] = particles[i].get_scalars();

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
	double diff = 1E300;	// Distance from peel_factor

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
		if(fabs(frac_below[i] - peel_factor) < diff)
		{
			which = i;
			diff = fabs(frac_below[i] - peel_factor);
		}
	}

	std::cout<<"# New threshold = ";
	for(size_t i=0; i<keep[which].size(); i++)
		std::cout<<keep[which][i]<<' ';
	std::cout<<std::endl;
	thresholds.push_back(keep[which]);

	// Write out dead points
	double log_dead_mass = log(frac_below[which]) + log_prior_mass;
	std::fstream fout("output.txt", std::ios::out|std::ios::app);
	fout<<std::setprecision(10);
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

	remove_redundant_thresholds();
}

