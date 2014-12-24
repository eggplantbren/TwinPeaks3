
/*********************************************************************
 *			IMPLEMENTATIONS BEGIN			     *
 *********************************************************************/

template<class Type>
Sampler<Type>::Sampler(int num_particles)
:num_particles(num_particles)
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
}

