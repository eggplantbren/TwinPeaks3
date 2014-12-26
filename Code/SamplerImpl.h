
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles)
:num_particles(num_particles)
,particles(num_particles)
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

	thresholds.clear();
}

template<class Type>
void Sampler<Type>::explore()
{
	const int steps = 10000;
	const int skip = 10;
	std::vector< std::vector<double> > keep(steps/skip);

	for(int i=0; i<steps; i++)
	{
		int which = DNest3::randInt(num_particles);
		Type proposal = particles[which];
		double logH = proposal.perturb();

		if(DNest3::randomU() <= exp(logH))
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
	double diff = 1E300;	// Distance from 1-exp(-1)

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
		if(fabs(frac_below[i] - (1. - exp(-1.))) < diff)
		{
			which = i;
			diff = fabs(frac_below[i] - (1. - exp(-1.)));
		}
	}

	std::cout<<"# New threshold = ";
	for(size_t i=0; i<keep[which].size(); i++)
		std::cout<<keep[which][i]<<' ';
	std::cout<<std::endl;
}

