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
,uccs(num_particles, std::vector<unsigned short>(num_particles))
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
    forbid_rectangles(worst, unique);
    replace_particle(worst);
    ++iteration;
}

template<class MyModel>
void Sampler<MyModel>::forbid_rectangles(size_t which, bool unique)
{
    // Shorter alias
    auto& rects = forbidden_rectangles;
    std::vector<ScalarType> latest(2);

    for(size_t i=0; i<num_particles; ++i)
    {
        for(int j=(num_particles-1); j>=0; --j)
        {
            if(uccs[i][j] == particle_uccs[which])
            {
//                std::cout<<"i = "<<i<<", j = "<<j<<std::endl;
                // The rectangle being added
                latest = {scalars[0][indices[0][j]], 
                          scalars[1][indices[1][num_particles-i-1]]};

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
                break;
            }
        }
    }


//    std::cout<<"Particles:"<<std::endl;
//    for(size_t i=0; i<num_particles; ++i)
//    {
//        std::cout<<std::setprecision(10);
//        std::cout<<scalars[0][i].get_value()<<' ';
//        std::cout<<scalars[1][i].get_value()<<' ';
//        std::cout<<particle_uccs[i]<<std::endl;
//    }

//    std::cout<<"Scalars:"<<std::endl;
//    for(size_t i=0; i<num_particles; ++i)
//    {
//        std::cout<<std::setprecision(10);
//        std::cout<<scalars[0][indices[1][i]].get_value()<<' ';
//        std::cout<<scalars[1][indices[1][i]].get_value()<<std::endl;
//    }


//    std::cout<<"Rectangles:"<<std::endl;
//    for(auto it=rects.begin(); it != rects.end(); ++it)
//    {
//        std::cout<<std::setprecision(10);
//        std::cout<<it->get_scalars()[0].get_value()<<' ';
//        std::cout<<it->get_scalars()[1].get_value()<<"     ";
//        std::cout<<it->get_opacity()<<std::endl;
//    }
//    exit(0);
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
    MyModel particle = particles[copy];
    std::vector<ScalarType> particle_scalars
                            {scalars[0][copy], scalars[1][copy]};
    double particle_ucc_tiebreaker = particle_ucc_tiebreakers[copy];
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
        proposal = particle;
        proposal_scalars = particle_scalars;
        proposal_ucc_tiebreaker = particle_ucc_tiebreaker;

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
            particle = proposal;
            particle_scalars = proposal_scalars;
            particle_ucc_tiebreaker = proposal_ucc_tiebreaker;
            logp = logp_proposal;
            ++accepted;
        }

//        std::cout<<logp<<' ';
//        for(int k=0; k<2; ++k)
//            std::cout<<particle_scalars[k].get_value()<<' ';
//        std::cout<<std::endl;
    }

    // Copy evolved particle back over to Sampler object arrays
    particles[which] = particle;
    scalars[0][which] = particle_scalars[0];
    scalars[1][which] = particle_scalars[1];
    particle_ucc_tiebreakers[which] = particle_ucc_tiebreaker;

    std::cout<<"# Done. Accepted "<<accepted<<"/"<<mcmc_steps<<" steps.";
    std::cout<<std::endl<<std::endl;
}

template<class MyModel>
double Sampler<MyModel>::log_prob(const std::vector<ScalarType>& s)
{
    // Alias
    const auto& rects = forbidden_rectangles;

    double logp = 0.0;
    for(auto it=rects.begin(); it != rects.end(); ++it)
    {
        if(it->get_opacity() >= 1.0)
        {
            if(ScalarType::compare(s, it->get_scalars()) == -1)
                return -std::numeric_limits<double>::max(); // -Infinity
        }
        else if(ScalarType::compare(s, it->get_scalars()) == -1)
            logp += log(1.0 - it->get_opacity());
    }
    return logp;
}

template<class MyModel>
void Sampler<MyModel>::calculate_uccs()
{
    // Start with zeros
    uccs.assign(num_particles,
                        std::vector<unsigned short>(num_particles, 0));

    // First construct the empirical measure. That is, put ones where particle
    // ranks are located, and zeroes elsewhere
    for(size_t i=0; i<num_particles; ++i)
        ++uccs[num_particles - ranks[1][i] - 1][ranks[0][i]];

//    for(size_t i=0; i<num_particles; ++i)
//    {
//        for(size_t j=0; j<num_particles; ++j)
//            std::cout<<uccs[i][j]<<' ';
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

    // Cumsum each row (increase to the left)
    for(size_t i=0; i<num_particles; ++i)
    {
        // Cumsum the row
        for(int j=(num_particles-2); j>=0; --j)
            uccs[i][j] += uccs[i][j+1];
    }
    // Cumsum each column (increase downwards)
    for(size_t j=0; j<num_particles; ++j)
    {
        // Cumsum the row
        for(size_t i=1; i<num_particles; ++i)
            uccs[i][j] += uccs[i-1][j];
    }

//    for(size_t i=0; i<num_particles; ++i)
//    {
//        for(size_t j=0; j<num_particles; ++j)
//            std::cout<<uccs[i][j]<<' ';
//        std::cout<<std::endl;
//    }
//    std::cout<<std::endl;

    // Assign the particle uccs
    for(size_t i=0; i<num_particles; ++i)
    {
        particle_uccs[i] = uccs[num_particles - ranks[1][i] - 1][ranks[0][i]];
//        std::cout<<particle_uccs[i]<<' ';
    }
//    std::cout<<std::endl;
}

} // namespace TwinPeaks

