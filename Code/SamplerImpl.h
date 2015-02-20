
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_threads, int num_particles, int steps, 					double peel_factor, int thin)
:num_threads(num_threads)
,num_particles(num_particles)
,steps(steps)
,thin(thin)
,peel_factor(peel_factor)
,rngs(num_threads)
,iterations(0)
,particles(num_particles)
,log_prior_mass(0.)
{
	// Initialise rngs
	for(int i=0; i<num_threads; i++)
		rngs[i] = gsl_rng_alloc(gsl_rng_mt19937);
	std::cout<<std::setprecision(8);
}

template<class Type>
int Sampler<Type>::badness(const Type& particle) const
{
	int count = 0;
	for(int i=int(thresholds.size())-1; i >= 0; i--)
	{
		if(is_below(particle.get_scalars(), thresholds[i],
				particle.get_tiebreakers(), thresholds_tiebreakers[i]))
		{
			count++;
			if(count == 2)
				return count;
		}
	}
	return count;
}

template<class Type>
Sampler<Type>::~Sampler()
{
	for(int i=0; i<num_threads; i++)
		gsl_rng_free(rngs[i]);
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
	thresholds_tiebreakers.clear();

	// Open output files
	std::fstream fout("output.txt", std::ios::out);
	fout.close();
	fout.open("sample.txt", std::ios::out);
	fout.close();

	// Set rng seeds
	for(int i=0; i<num_threads; i++)
		gsl_rng_set(rngs[i], time(0) + 10*(i+1));
}

template<class Type>
void Sampler<Type>::runThread(int thread, const std::vector<int>& which_particles, std::vector<int>& bad, int& accepts)
{
	// Initialise thread-specific rng
	DNest3::RandomNumberGenerator::initialise_instance();
	DNest3::RandomNumberGenerator::get_instance().set_rng(rngs[thread]);

	int which, proposal_badness;
	Type proposal;
	double logH;

	for(int i=0; i<steps; i++)
	{
		which = which_particles[DNest3::randInt(which_particles.size())];

		proposal = particles[which];
		logH = proposal.perturb();
		proposal.perturb_tiebreakers();
		proposal_badness = badness(proposal);

		if(proposal_badness <= bad[which] &&
				DNest3::randomU() <= exp(logH))
		{
			particles[which] = proposal;
			bad[which] = proposal_badness;
			accepts++;
		}
	}

	// Save thread-specific rng
	gsl_rng_memcpy(rngs[thread], DNest3::RandomNumberGenerator::get_instance().get_rng());
}

// Exploration just to ensure everything's above all thresholds
template<class Type>
void Sampler<Type>::refresh()
{
	std::cout<<"# Evolving particles..."<<std::flush;

	std::vector<int> bad(num_particles);
	std::vector<int> need_refresh;
	for(int i=0; i<num_particles; i++)
	{
		bad[i] = badness(particles[i]);
		if(bad[i] != 0)
			need_refresh.push_back(i);
	}

	// Divide up the work between the threads
	std::vector< std::vector<int> > which_particles;
	for(size_t i=0; i<need_refresh.size(); i++)
	{
		if(int(which_particles.size()) < num_threads)
			which_particles.push_back(std::vector<int>());
		which_particles[i%num_threads].push_back(need_refresh[i]);
	}
	std::vector<int> accepts(which_particles.size(), 0);

//	// Create the threads
//	boost::thread_group threads;
//	for(size_t i=0; i<which_particles.size(); i++)
//	{
//		threads.create_thread(boost::bind(&Sampler<Type>::runThread, this, i, boost::ref(which_particles[i]), boost::ref(bad), boost::ref(accepts[i])));
//	}
//	threads.join_all();

	// Make sure they're not all bad
	bool all_bad = true;
	for(int i=0; i<num_particles; i++)
		if(bad[i] == 0)
			all_bad = false;

	if(all_bad)
		std::cerr<<"# WARNING: All particles are bad."<<std::endl;

	// Resample any bad points by copying good ones
	int copy;
	int num_copied = 0;
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
			num_copied++;
		}
	}

	// Create the threads
	boost::thread_group threads2;
	for(size_t i=0; i<which_particles.size(); i++)
	{
		threads2.create_thread(boost::bind(&Sampler<Type>::runThread, this, i, boost::ref(which_particles[i]), boost::ref(bad), boost::ref(accepts[i])));
	}
	threads2.join_all();

	int total_accepts = 0;
	for(size_t i=0; i<accepts.size(); i++)
		total_accepts += accepts[i];

	std::cout<<"done. Acceptance rate = "<<total_accepts<<"/"<<(steps*accepts.size())<<std::endl<<"."<<std::endl;
	remove_redundant_thresholds();
}

template<class Type>
void Sampler<Type>::remove_redundant_thresholds()
{
	top:
	for(size_t i=0; i<thresholds.size()-1; i++)
	{
		if(is_below(thresholds[i], thresholds.back(), thresholds_tiebreakers[i], thresholds_tiebreakers.back()))
		{
			thresholds.erase(thresholds.begin() + i);
			thresholds_tiebreakers.erase(thresholds_tiebreakers.begin() + i);
			goto top;
		}
	}
}

template<class Type>
void Sampler<Type>::explore()
{
	std::vector< std::vector<double> > keep(num_particles);
	std::vector< std::vector<double> > keep_tiebreakers(num_particles);

	for(int i=0; i<num_particles; i++)
	{
		keep[i] = particles[i].get_scalars();
		keep_tiebreakers[i] = particles[i].get_tiebreakers();
	}

	create_threshold(keep, keep_tiebreakers);
}

template<class Type>
bool Sampler<Type>::is_below(const std::vector<double>& s1,
				const std::vector<double>& s2,
				const std::vector<double>& tb1,
				const std::vector<double>& tb2) const
{
	for(size_t i=0; i<s1.size(); i++)
	{
		if(s1[i] > s2[i] || (s1[i] == s2[i] && tb1[i] >= tb2[i]))
			return false;
	}
	return true;
}

template<class Type>
void Sampler<Type>::create_threshold(const std::vector< std::vector<double> >&
						keep,
					const std::vector< std::vector<double> >&
						keep_tiebreakers)
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
				frac_below[i] += is_below(keep[j], keep[i],
							keep_tiebreakers[j],
							keep_tiebreakers[i]);
		}
		frac_below[i] /= (keep.size() - 1);
		if(fabs(log(frac_below[i]+1E-300) - log(peel_factor)) < diff)
		{
			which = i;
			diff = fabs(log(frac_below[i]+1E-300) - log(peel_factor));
		}
	}
	// All candidates to become the next rectangle
	// Choose uniformly from them
	std::vector<int> indices;
	for(size_t i=0; i<keep.size(); i++)
		if(fabs(log(frac_below[i]+1E-300) - log(peel_factor)) == diff)
			indices.push_back(i);
	which = indices[DNest3::randInt(indices.size())];

	std::cout<<"# New threshold = ";
	for(size_t i=0; i<keep[which].size(); i++)
		std::cout<<keep[which][i]<<' ';
	std::cout<<std::endl;
	thresholds.push_back(keep[which]);
	thresholds_tiebreakers.push_back(keep_tiebreakers[which]);

	// Write out dead points
	double log_dead_mass = log(frac_below[which]) + log_prior_mass;
	std::fstream fout("output.txt", std::ios::out|std::ios::app);
	fout<<std::setprecision(10);
	std::vector<int> dead;
	for(size_t i=0; i<keep.size(); i++)
	{
		if(int(i) != which && is_below(keep[i], keep[which], keep_tiebreakers[i], keep_tiebreakers[which]))
		{
			fout<<(log_dead_mass - log(keep.size()*frac_below[which]))<<' ';
			for(size_t j=0; j<keep[i].size(); j++)
				fout<<keep[i][j]<<' ';
			fout<<std::endl;
			dead.push_back(i);
		}
	}
	fout.close();

	log_prior_mass = DNest3::logdiffexp(log_prior_mass, log_dead_mass);
	std::cout<<"# Peeling away "<<frac_below[which]<<" of the remaining prior mass."<<std::endl;
	std::cout<<"# log(remaining prior mass) = "<<log_prior_mass<<std::endl;

	iterations++;
	// Search for a bad particle to write to file, along with an importance
	// weight
	if(iterations%thin == 0)
	{
		fout.open("sample.txt", std::ios::out|std::ios::app);
		if(dead.size() == 0)
			std::cerr<<"# WARNING: Couldn't find a dead particle to write to disk."<<std::endl;
		else
		{
			int which = DNest3::randInt(dead.size());
			fout<<log_dead_mass<<' ';
			fout<<particles[dead[which]]<<std::endl;
		}
		fout.close();
	}

//	remove_redundant_thresholds();
}

