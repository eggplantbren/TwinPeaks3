#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <limits>
#include <fstream>
#include <algorithm>
#include <cassert>
#include <thread>
#include <set>
#include "Utils.h"

namespace TwinPeaks
{

template<class MyModel>
Sampler<MyModel>::Sampler(const std::vector<RNG>& rngs, int num_particles,
                            int mcmc_steps,    int saves_per_iteration)
:rngs(rngs)
,num_particles(num_particles)
,particles(num_particles)
,scalar1(num_particles)
,scalar2(num_particles)
,scalar1_sorted(num_particles)
,scalar2_sorted(num_particles)
,uccs(num_particles, std::vector<unsigned short>(num_particles))
,particle_uccs(num_particles)
,status(num_particles)
,mcmc_steps(mcmc_steps)
,saves_per_iteration(saves_per_iteration)
,initialised(false)
,iteration(0)
,log_prior_mass(0.0)
{
    // Open and close outputs file to clear them
    std::fstream fout;
    fout.open("sample.txt", std::ios::out);
    fout.close();
    fout.open("sample_info.txt", std::ios::out);
    fout.close();
}

template<class MyModel>
void Sampler<MyModel>::initialise()
{
    std::cout<<"# Generating "<<num_particles<<" particles from the prior...";
    std::cout<<std::flush;
    for(int i=0; i<num_particles; ++i)
    {
        particles[i].from_prior(rngs[0]);
        auto& s = particles[i].get_scalars();
        scalar1[i] = s[0];
        scalar2[i] = s[1];
        scalar1[i].from_prior(rngs[0]);
        scalar2[i].from_prior(rngs[0]);
    }
    std::cout<<"done."<<std::endl<<std::endl;
}

template<class MyModel>
void Sampler<MyModel>::calculate_uccs()
{
    // Agsort the two scalars
    std::vector<size_t> indices1 = argsort(scalar1);
    std::vector<size_t> indices2 = argsort(scalar2);

    // Do the sorting
    for(int i=0; i<num_particles; i++)
    {
        scalar1_sorted[i] = scalar1[indices1[i]];
        scalar2_sorted[i] = scalar2[indices2[i]];
    }

    // Calculate ranks
    std::vector<size_t> r1 = compute_ranks(indices1);
    std::vector<size_t> r2 = compute_ranks(indices2);

    // Calculate empirical measure of ranks
    std::vector< std::vector<bool> > empirical_measure(num_particles,
                                        std::vector<bool>(num_particles, 0));
    for(int i=0; i<num_particles; i++)
        empirical_measure[num_particles - r2[i] - 1][r1[i]] = 1;

    // Integrate empirical measure to get (inclusive) upper corner count
    uccs.assign(num_particles, std::vector<unsigned short>(num_particles, 0));
    int up, right, up_and_right;
    for(int i=0; i<num_particles; ++i)
    {
        for(int j=(num_particles-1); j>=0; --j)
        {
            up = 0;
            right = 0;
            up_and_right = 0;
            if(i != 0)
                up = uccs[i-1][j];
            if(j != (num_particles-1))
                right = uccs[i][j+1];
            if((i != 0) && (j != (num_particles-1)))
                up_and_right = uccs[i-1][j+1];
            uccs[i][j] = empirical_measure[i][j] + up + right - up_and_right;
        }
    }

    // The uccs of the particles themselves
    for(int i=0; i<num_particles; ++i)
        particle_uccs[i] = uccs[num_particles - r2[i] - 1][r1[i]];
}

template<class MyModel>
unsigned short Sampler<MyModel>::choose_ucc_threshold()
{
    // Sort the particle uccs from highest to lowest
    auto particle_uccs_sorted = particle_uccs;
    std::sort(particle_uccs_sorted.begin(), particle_uccs_sorted.end());
    std::reverse(particle_uccs_sorted.begin(), particle_uccs_sorted.end());

    // Ensure there are at least three unique ucc values
    short unique = 1;
    for(int i=1; i<num_particles; ++i)
        if(particle_uccs_sorted[i] != particle_uccs_sorted[i-1])
            ++unique;
    if(unique < 3)
    {
        std::cerr<<"# Less than three distinct particle_ucc values."<<std::endl;
        exit(0);
    }

    // Find a ucc threshold
    std::vector<int> threshold_id_candidates;
    for(int i=1; i<(num_particles-1); ++i)
        if(particle_uccs_sorted[i] != particle_uccs_sorted[i-1])
            threshold_id_candidates.push_back(i);

    // Find the threshold that's closest to dividing the particles in half
    std::vector<double> balance(threshold_id_candidates.size());
    unsigned int threshold;
    for(size_t i=0; i<threshold_id_candidates.size(); ++i)
    {
        threshold = particle_uccs_sorted[threshold_id_candidates[i]];

        // Count interior and exterior particles
        int num_interior = 0;
        int num_exterior = 0;
        for(int j=0; j<num_particles; ++j)
        {
            if(particle_uccs[j] < threshold)
                ++num_exterior;
            if(particle_uccs[j] > threshold)
                ++num_interior;
        }
        balance[i] = static_cast<double>(num_interior + 1E-12)/
                                        (num_interior + num_exterior + 1E-12);
        balance[i] = log(balance[i]/(1.0 - balance[i]));
    }
    int best = 0;
    for(size_t i=1; i<balance.size(); ++i)
        if(std::abs(balance[i]) < std::abs(balance[best]))
            best = i;
    threshold = particle_uccs_sorted[threshold_id_candidates[best]];

    // Calculate particle statuses
    for(int j=0; j<num_particles; ++j)
    {
        if(particle_uccs[j] < threshold)
            status[j] = 1;
        if(particle_uccs[j] == threshold)
            status[j] = 0;
        if(particle_uccs[j] > threshold)
            status[j] = -1;
    }

    return threshold;
}

template<class MyModel>
double Sampler<MyModel>::do_iteration()
{
    // Calculate the uccs and choose a threshold
    calculate_uccs();
    auto threshold = choose_ucc_threshold();

    // Place forbidding rectangles anywhere ucc >= threshold
    std::vector< std::vector<ScalarType> > latest_rectangles;
    for(int i=0; i<num_particles; ++i)
    {
        int j = num_particles-1;
        while(j > 0 && uccs[i][j] < threshold)
            --j;
        if(uccs[i][j] >= threshold)
        {
            std::vector<ScalarType> latest{scalar1_sorted[j],
                                           scalar2_sorted[num_particles-i-1]};
            latest_rectangles.push_back(latest);
        }
    }
    context.add_opaque_rectangles(latest_rectangles);

    // Count number of particles with each status
    int num_interior = 0;
    int num_boundary = 0;
    int num_exterior = 0;
    for(int i=0; i<num_particles; ++i)
    {
        if(status[i] == -1)
            ++num_interior;
        else if(status[i] == 0)
            ++num_boundary;
        else
            ++num_exterior;
    }

    // Save interior particles
    for(int i=0; i<num_particles; ++i)
    {
        if(status[i] == -1)
        {
            // Assign prior weight
            double logw = log_prior_mass - log(num_interior + num_exterior);

            // Write it out to an output file
            std::fstream fout;
            fout.open("sample.txt", std::ios::out|std::ios::app);
            fout<<logw<<' ';
            fout<<scalar1[i].get_value()<<' '<<scalar2[i].get_value()<<' ';
            particles[i].write_text(fout);
            fout<<std::endl;
            fout.close();

            fout.open("sample_info.txt", std::ios::out|std::ios::app);
            fout<<logw<<' ';
            fout<<scalar1[i].get_value()<<' '<<scalar2[i].get_value()<<' ';
            fout<<std::endl;
            fout.close();
        }
    }

    // Reduce remaining prior mass
    double surviving_fraction = static_cast<double>(num_exterior)/
                                                (num_interior + num_exterior);
    log_prior_mass += log(surviving_fraction);

    // Print messages
    std::cout<<"# Iteration "<<(iteration+1)<<"."<<std::endl;
    std::cout<<"# (num_interior, num_boundary, num_exterior) = (";
    std::cout<<num_interior<<", "<<num_boundary<<", "<<num_exterior<<").";
    std::cout<<std::endl;
    std::cout<<"# Killing "<<(num_interior + num_boundary)<<" particles. "<<std::endl;
    std::cout<<"# There are now "<<context.get_num_rectangles();
    std::cout<<" rectangles, and log(remaining prior mass) = ";
    std::cout<<log_prior_mass<<"."<<std::endl;

    replace_dead_particles();

    ++iteration;
    return log_prior_mass;
}


template<class MyModel>
void Sampler<MyModel>::replace_dead_particles()
{
    // Replace dead particles
    int num_threads = rngs.size();

	// Backup
    backup_particles = particles;
    backup_scalar1 = scalar1;
    backup_scalar2 = scalar2;

    // Assign particles to threads
    std::vector< std::vector<int> > which_particles(num_threads);
    int k=0;
    int count = 0;
    for(int i=0; i<num_particles; ++i)
    {
        if(status[i] != 1)
        {
            ++count;
            which_particles[(k++)%num_threads].push_back(i);
        }
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

    std::cout<<"# Accepted "<<accepted<<"/"<<count*mcmc_steps<<" (";
    std::cout<<std::fixed<<std::setprecision(1);
    std::cout<<(100.*accepted/(count*mcmc_steps));
    std::cout<<"%)."<<std::endl<<std::endl;
    std::cout<<std::defaultfloat<<std::setprecision(6);
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
    scalar1[which] = backup_scalar1[copy];
    scalar2[which] = backup_scalar2[copy];

    // Do the MCMC
    MyModel proposal;
    std::vector<ScalarType> s_proposal;
    double logH;
    int accepted = 0;
    for(int i=0; i<mcmc_steps; i++)
    {
        proposal = particles[which];
        s_proposal = {scalar1[which], scalar2[which]};

        logH = proposal.perturb(rngs[which_rng]);
        for(size_t j=0; j<s_proposal.size(); j++)
            s_proposal[j].set_value(proposal.get_scalars()[j]);
        logH += s_proposal[rngs[which_rng].rand_int(s_proposal.size())].perturb(rngs[which_rng]);

        if(rngs[which_rng].rand() <= exp(logH) &&
           context.log_prob(s_proposal) == 0.0)
        {
            particles[which] = proposal;
            scalar1[which] = s_proposal[0];
            scalar2[which] = s_proposal[1];
            ++accepted;
        }
    }
    return accepted;
}

} // namespace TwinPeaks


