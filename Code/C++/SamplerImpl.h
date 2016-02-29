#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include "Utils.h"

namespace TwinPeaks
{

template<class MyModel>
Sampler<MyModel>::Sampler(unsigned int num_particles, unsigned int mcmc_steps,
                            const std::vector<RNG>& rngs)
:num_particles(num_particles)
,mcmc_steps(mcmc_steps)
,iteration(1)
,rngs(rngs)
,particles(num_particles)
,scalars(2, std::vector<ScalarType>(num_particles))
,indices(2, std::vector<size_t>(num_particles))
,ranks(2, std::vector<size_t>(num_particles))
,particle_uccs(num_particles)
,particle_ucc_tiebreakers(num_particles)
{
    assert(num_particles > 0 && mcmc_steps > 0);
}

template<class MyModel>
void Sampler<MyModel>::initialise()
{
    std::cout<<"# Initialising TwinPeaks sampler."<<std::endl;
    std::cout<<"# Generating "<<num_particles<<" particles from the prior...";
    std::cout<<std::flush;
    for(size_t k=0; k<num_particles; ++k)
    {
        particles[k].from_prior(rngs[0]);
        const std::vector<double>& s = particles[k].get_scalars();
        for(size_t i=0; i<2; ++i)
            scalars[i][k] = ScalarType(s[i], rngs[0].rand());
        particle_ucc_tiebreakers[k] = rngs[0].rand();
    }
    std::cout<<"done."<<std::endl<<std::endl;
}

template<class MyModel>
void Sampler<MyModel>::do_iteration()
{
    // Argsort by the two scalars
    for(int i=0; i<2; ++i)
    {
        indices[i] = argsort(scalars[i]);
        ranks[i] = compute_ranks(indices[i]);
    }

    // Calculate the UCCs
    calculate_uccs();

    // Find worst ucc (highest!)
    size_t worst = 0;
    for(size_t i=1; i<num_particles; ++i)
        if(particle_uccs[i] > particle_uccs[worst] ||
            (
            particle_uccs[i] == particle_uccs[worst] &&
            particle_ucc_tiebreakers[i] > particle_ucc_tiebreakers[worst]))
            worst = i;

    // Check to see if the worst particle ucc is unique
    bool unique = true;
    for(size_t i=0; i<num_particles; ++i)
    {
        if(i != worst && particle_uccs[i] == particle_uccs[worst])
        {
            unique = false;
            break;
        }
    }

    // Compression
    double log_X = -static_cast<double>(iteration)/num_particles;

    // Print information to the screen
    std::cout<<"# Iteration "<<iteration<<". log(X) = "<<log_X<<std::endl;

    // Open the output file
    std::fstream fout;
    if(iteration == 1)
        fout.open("output.txt", std::ios::out);
    else
        fout.open("output.txt", std::ios::out | std::ios::app);

    // Write info to disk
    fout<<std::setprecision(12);
    fout<<log_X<<' ';
    fout<<scalars[0][worst].get_value()<<' ';
    fout<<scalars[1][worst].get_value()<<std::endl;

    // Close output file
    fout.close();

    // Restrict the space and generate a replacement particle
    forbid_rectangle(worst, unique);
    replace_particle(worst);
    ++iteration;
}

template<class MyModel>
void Sampler<MyModel>::forbid_rectangle(size_t which, bool unique)
{
    // Shorter alias
    auto& rects = forbidden_rectangles;

    // The rectangle being added
    std::vector<ScalarType> latest{scalars[0][which], scalars[1][which]};

    // Remove redundant rectangles
    if(unique)
    {
        for(auto it=rects.begin(); it != rects.end(); ++it)
        {
	        if(ScalarType::compare(latest, it->get_scalars()) == 1)
		        it = rects.erase(it);
        }
    }

    // Forbid the rectangle
    if(unique)
        rects.push_front(Rectangle(latest, 1.0));
    else
        rects.push_front(Rectangle(latest, particle_ucc_tiebreakers[which]));
}

template<class MyModel>
void Sampler<MyModel>::replace_particle(size_t which)
{
    std::cout<<"# There are now "<<forbidden_rectangles.size();
    std::cout<<" forbidden rectangles."<<std::endl;
    std::cout<<"# Replacing particle by cloning and doing MCMC..."<<std::endl;

    // Choose another particle to clone
    size_t copy;
    do
    {
        copy = rngs[0].rand_int(num_particles);
    }while(copy == which);

    // Clone it
    particles[which] = particles[copy];
    std::vector<ScalarType> particle_scalars
                            {scalars[0][copy], scalars[1][copy]};
    double logp = log_prob(particle_scalars);

    // Variables that we'll need to propose
    MyModel proposal;
    std::vector<ScalarType> proposal_scalars(2);
    double proposal_ucc_tiebreaker;
    double logp_proposal;

    // Counter for M-H acceptance fraction
    size_t accepted = 0;

    // Do MCMC
    double log_H;
    for(size_t i=0; i<mcmc_steps; ++i)
    {
        // Copy to do proposal
        proposal = particles[which];
        proposal_scalars = particle_scalars;
        proposal_ucc_tiebreaker = particle_ucc_tiebreakers[which];

        // Propose
        log_H = proposal.perturb(rngs[0]);
        for(size_t j=0; j<2; ++j)
        {
            proposal_scalars[j].set_value(proposal.get_scalars()[j]);
            proposal_scalars[j].perturb(rngs[0]);
        }
        proposal_ucc_tiebreaker += rngs[0].randh();
        wrap(proposal_ucc_tiebreaker, 0.0, 1.0);

        // Measure the target density for the proposed point
        logp_proposal = log_prob(proposal_scalars);

        // Accept
        if(rngs[0].rand() <= exp(logp_proposal - logp + log_H))
        {
            particles[which] = proposal;
            particle_scalars = proposal_scalars;
            particle_ucc_tiebreakers[which] = proposal_ucc_tiebreaker;
            ++accepted;
        }
    }

    scalars[0][which] = particle_scalars[0];
    scalars[1][which] = particle_scalars[1];

    std::cout<<"# Done. Accepted "<<accepted<<"/"<<mcmc_steps<<" steps.";
    std::cout<<std::endl<<std::endl;
}

template<class MyModel>
double Sampler<MyModel>::log_prob(const std::vector<ScalarType>& s)
{
    double logp = 0.0;
    for(const auto& rect: forbidden_rectangles)
    {
        if(rect.get_opacity() >= 1.0)
        {
            if(ScalarType::compare(s, rect.get_scalars()) == -1)
                return -std::numeric_limits<double>::max(); // -Infinity
        }
        else if(ScalarType::compare(s, rect.get_scalars()) == -1)
            logp += log(1.0 - rect.get_opacity());
    }
    return logp;
}

template<class MyModel>
void Sampler<MyModel>::calculate_uccs()
{
    // Start with zeros
    std::vector< std::vector<unsigned short> > uccs(num_particles,
                        std::vector<unsigned short>(num_particles, 0));

    // First construct the empirical measure. That is, put ones where particle
    // ranks are located, and zeroes elsewhere
    for(size_t i=0; i<num_particles; ++i)
        ++uccs[ranks[0][i]][ranks[1][i]];

    // Sum over horizontal direction
    for(size_t i=0; i<num_particles; ++i)
    {
        for(int j=num_particles-2; j>=0; --j)
            uccs[i][j] += uccs[i][j+1];
    }
    // Sum over vertical direction
    for(size_t i=1; i<num_particles; ++i)
    {
        for(size_t j=0; j<num_particles; ++j)
            uccs[i][j] += uccs[i-1][j];
    }

    // Assign the particle uccs (subtract 1 to not count self)
    for(size_t i=0; i<num_particles; ++i)
        particle_uccs[i] = uccs[ranks[0][i]][ranks[1][i]] - 1;
}

/*
template<class MyModel>
double Sampler<MyModel>::do_iteration()
{
	// Extract and sort the two scalars
	std::vector<ScalarType> s1(num_particles);
	std::vector<ScalarType> s2(num_particles);
	for(int i=0; i<num_particles; i++)
	{
		s1[i] = scalars[i][0];
		s2[i] = scalars[i][1];
	}
	// Calculate ranks
	std::vector<size_t> r1 = ranks(s1);
	std::vector<size_t> r2 = ranks(s2);
	// Sort
	std::sort(s1.begin(), s1.end());
	std::sort(s2.begin(), s2.end());

	// Calculate empirical measure of ranks
	std::vector< std::vector<bool> > empirical_measure(num_particles,
										std::vector<bool>(num_particles, 0));
	for(int i=0; i<num_particles; i++)
		empirical_measure[num_particles - r2[i] - 1][r1[i]] = 1;

	// Integrate empirical measure to get (inclusive) upper corner count
	std::vector< std::vector<unsigned short> > ucc(num_particles,
										std::vector<unsigned short>(num_particles, 0));
	int up, right, up_and_right;
	for(int i=0; i<num_particles; i++)
	{
		for(int j=(num_particles-1); j>=0; j--)
		{
			up = 0;
			right = 0;
			up_and_right = 0;
			if(i != 0)
				up = ucc[i-1][j];
			if(j != (num_particles-1))
				right = ucc[i][j+1];
			if((i != 0) && (j != (num_particles-1)))
				up_and_right = ucc[i-1][j+1];
			ucc[i][j] = empirical_measure[i][j] + up + right - up_and_right;
		}
	}

	// The uccs of the particles themselves
	std::vector<int> particle_uccs(num_particles);
	for(int i=0; i<num_particles; i++)
		particle_uccs[i] = ucc[num_particles - r2[i] - 1][r1[i]];

	// Sort the particle uccs from highest to lowest
	std::vector<int> particle_uccs_sorted = particle_uccs;
	std::sort(particle_uccs_sorted.begin(), particle_uccs_sorted.end());
	std::reverse(particle_uccs_sorted.begin(), particle_uccs_sorted.end());

	// Make a ucc threshold (particles on the threshold die too)
	int threshold, threshold_id;
	int count_thresholds_tried = 0;

	threshold_id = num_particles/2;
	threshold_selection:
	threshold = particle_uccs_sorted[threshold_id];
	++count_thresholds_tried;

	// -1 <===> interior
	//  0 <===> boundary
	// +1 <===> exterior
	status.assign(num_particles, 0);
	int num_interior = 0;
	int num_exterior = 0;
	int num_boundary = 0;
	for(int i=0; i<num_particles; i++)
	{
		if(particle_uccs[i] > threshold)
		{
			status[i] = -1;
			++num_interior;
		}
		if(particle_uccs[i] < threshold)
		{
			status[i] = +1;
			++num_exterior;
		}
	}
	num_boundary = num_particles - (num_interior + num_exterior);

	// Handle the case where there are no interior or exterior particles
	if(num_interior == 0 || num_exterior == 0)
	{
		if(count_thresholds_tried == num_particles)
		{
			// Print messages and quit
			std::cout<<"# Iteration "<<(iteration+1)<<"."<<std::endl;
			std::cout<<"# (num_interior, num_boundary, num_exterior) = (";
			std::cout<<num_interior<<", "<<num_boundary<<", "<<num_exterior<<")."<<std::endl;
			std::cout<<"# CANNOT CONTINUE."<<std::endl;
			exit(0);
		}

		if(threshold_id < num_particles/2 && threshold_id != (num_particles-1))
			--threshold_id;
		else if(threshold_id == (num_particles-1))
			threshold_id = num_particles/2 - 1;
		else
			++threshold_id;

		goto threshold_selection;
	}
	else // Standard TwinPeaks
	{
		// Place forbidding rectangles anywhere ucc >= threshold
		for(int i=0; i<num_particles; i++)
		{
			int j = num_particles-1;
			while(j > 0 && ucc[i][j] < threshold)
				j--;

			if(ucc[i][j] >= threshold)
			{
				std::vector<ScalarType> latest{s1[j], s2[num_particles-i-1]};
				prune_rectangles(latest);
				rects.push_front(latest);
			}
		}
	}

	// Select some particles to save in their entirety
	std::vector<bool> save(num_particles, false);
	if(num_interior != 0)
	{
		for(int i=0; i<saves_per_iteration; i++)
		{
			int ii;
			do
			{
				ii = rngs[0].rand_int(num_particles);
			}while(status[ii] != -1);
			save[ii] = true;
		}
	}

	// Save interior particles
	for(int i=0; i<num_particles; i++)
	{
		if(status[i] == -1)
		{
			// Assign prior weight
			double logw = log_prior_mass - log(num_particles - num_boundary);

			// Write it out to an output file
			std::fstream fout;
			if(save[i])
			{
				fout.open("sample.txt", std::ios::out|std::ios::app);
				fout<<logw<<' ';
				for(ScalarType s: scalars[i])
					fout<<s.get_value()<<' ';
				particles[i].write_text(fout);
				fout<<std::endl;
				fout.close();
			}
			fout.open("sample_info.txt", std::ios::out|std::ios::app);
			fout<<logw<<' ';
			for(ScalarType s: scalars[i])
				fout<<s.get_value()<<' ';
			fout<<std::endl;
			fout.close();
		}
	}

	// Reduce remaining prior mass
	log_prior_mass = logdiffexp(log_prior_mass, log_prior_mass + log(num_interior) - log(num_particles - num_boundary));

	// Print messages
	std::cout<<"# Iteration "<<(iteration+1)<<"."<<std::endl;
	std::cout<<"# (num_interior, num_boundary, num_exterior) = (";
	std::cout<<num_interior<<", "<<num_boundary<<", "<<num_exterior<<")."<<std::endl;
	std::cout<<"# Killing "<<(num_interior + num_boundary)<<" particles. "<<std::endl;
	std::cout<<"# "<<rects.size()<<" rectangles. Log(remaining prior mass) = ";
	std::cout<<log_prior_mass<<"."<<std::endl;

	// Replace dead particles
	int num_threads = rngs.size();

	// Backup
	backup_particles = particles;
	backup_scalars = scalars;

	// Assign particles to threads
	std::vector< std::vector<int> > which_particles(num_threads);
	int k=0;
	for(int i=0; i<num_particles; ++i)
	{
		if(status[i] != 1)
			which_particles[(k++)%num_threads].push_back(i);
	}
	// Store acceptance counts
	std::vector<int> accepts(num_threads);

	// Do MCMC to generate new particles
	std::vector<std::thread> threads;
	for(int i=0; i<num_threads; i++)
	{
		threads.emplace_back(std::thread(std::bind(&Sampler<MyModel>::refresh_particles, this, which_particles[i], i, std::ref(accepts[i]))));
	}
	for(int i=0; i<num_threads; i++)
		threads[i].join();
	// Sum acceptance counts
	int accepted = 0;
	for(const int& c: accepts)
		accepted += c;

	std::cout<<"# Accepted "<<accepted<<"/"<<(num_interior + num_boundary)*mcmc_steps<<" (";
	std::cout<<std::fixed<<std::setprecision(1);
	std::cout<<(100.*accepted/((num_interior + num_boundary)*mcmc_steps));
	std::cout<<"%)."<<std::endl<<std::endl;
	std::cout<<std::defaultfloat<<std::setprecision(6);

	iteration++;
	return log_prior_mass;
}


// Refresh all the particles in indices
template<class MyModel>
void Sampler<MyModel>::refresh_particles(const std::vector<int>& indices,
														int which_rng,
														int& accepts)
{
	accepts = 0;
	for(int i: indices)
		accepts += refresh_particle(i, which_rng);
}

template<class MyModel>
int Sampler<MyModel>::refresh_particle(int which, int which_rng)
{
	// Choose a particle to clone
	int copy;
	do
	{
		copy = rngs[which_rng].rand_int(num_particles);
	}
	while(status[copy] != 1);

	// Clone it
	particles[which] = backup_particles[copy];
	scalars[which] = backup_scalars[copy];

	// Do the MCMC
	MyModel proposal;
	std::vector<ScalarType> s_proposal;
	double logH;
	int accepted = 0;
	for(int i=0; i<mcmc_steps; i++)
	{
		proposal = particles[which];
		s_proposal = scalars[which];

		logH = proposal.perturb(rngs[which_rng]);
		for(size_t j=0; j<s_proposal.size(); j++)
			s_proposal[j].set_value(proposal.get_scalars()[j]);
		logH += s_proposal[rngs[which_rng].rand_int(s_proposal.size())].perturb(rngs[which_rng]);

		if(rngs[which_rng].rand() <= exp(logH) && is_okay(s_proposal))
		{
			particles[which] = proposal;
			scalars[which] = s_proposal;
			accepted++;
		}
	}
	return accepted;
}
*/

} // namespace TwinPeaks

