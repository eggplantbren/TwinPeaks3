#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <limits>
#include "Utils.h"

template<class MyModel>
Results<MyModel>::Results()
:read_sample_info_flag(false)
{

}

template<class MyModel>
void Results<MyModel>::read_sample_info()
{
	// Open the output files (assume standard filenames)
	std::fstream sample_info_file("sample_info.dat", std::ios::in|std::ios::binary);
	if(!sample_info_file)
	{
		std::cerr<<"Error opening sample_info.dat.";
		std::cerr<<std::endl;
		return;
	}

	// Clear the values in memory
	logX.clear();
	log_prior_masses.clear();
	log_likelihoods.clear();

	// Read in the number of particles
	sample_info_file.read(reinterpret_cast<char*>(&num_particles),
									sizeof(num_particles));

	// Read until end of file
	int temp1; double temp2, temp3, temp4;
	while(sample_info_file.read(reinterpret_cast<char*>(&temp1), sizeof(temp1)))
	{
		sample_info_file.read(reinterpret_cast<char*>(&temp2), sizeof(temp2));
		sample_info_file.read(reinterpret_cast<char*>(&temp3), sizeof(temp3));
		sample_info_file.read(reinterpret_cast<char*>(&temp4), sizeof(temp4));

		// Keep the prior and (unnormalised) posterior masses in memory
		logX.push_back(temp2);
		log_likelihoods.push_back(temp3);
	}
	sample_info_file.close();

	// Calculate prior masses
	calculate_prior_masses();

	// Calculate log evidence and information
	calculate_log_evidence();
	calculate_information();

	read_sample_info_flag = true;
}

template<class MyModel>
void Results<MyModel>::calculate_prior_masses()
{
	log_prior_masses.resize(logX.size());
	log_prior_masses[0] = logdiffexp(0., logX[0]);
	for(size_t i=1; i<log_prior_masses.size(); i++)
		log_prior_masses[i] = logdiffexp(logX[i-1], logX[i]);
}

template<class MyModel>
void Results<MyModel>::calculate_log_evidence()
{
	logZ = -std::numeric_limits<double>::max();
	for(size_t i=0; i<log_prior_masses.size(); i++)
		logZ = logsumexp(logZ, log_prior_masses[i] + log_likelihoods[i]);
}

template<class MyModel>
void Results<MyModel>::calculate_information()
{
	H = 0.;

	double logp;
	for(size_t i=0; i<log_prior_masses.size(); i++)
	{
		// Normalised posterior mass
		logp = log_prior_masses[i] + log_likelihoods[i] - logZ;
		H += exp(logp)*(logp - log_prior_masses[i]);
	}
}

template<class MyModel>
void Results<MyModel>::regenerate_logX(RNG& rng)
{
	// Using the method of generating rands and taking their max
	std::vector<double> logx(num_particles);
	for(size_t i=0; i<logx.size(); i++)
		logx[i] = log(rng.rand());

	std::vector<double>::iterator winner;
	for(size_t i=0; i<logX.size(); i++)
	{
		winner = max_element(logx.begin(), logx.end());
		logX[i] = *winner;
		*winner += log(rng.rand());
	}

	calculate_prior_masses();
	calculate_log_evidence();
	calculate_information();
}

template<class MyModel>
void Results<MyModel>::generate_posterior_samples() const
{
	std::cout<<"# Generating posterior samples. Effective sample size = ";

	// Calculate effective sample size
	double ESS = 0.;

	double logp;
	for(size_t i=0; i<log_prior_masses.size(); i++)
	{
		// Normalised posterior mass
		logp = log_prior_masses[i] + log_likelihoods[i] - logZ;
		ESS += -exp(logp)*logp;
	}
	ESS = exp(ESS);
	std::cout<<ESS<<"."<<std::endl;

	// Open sample.dat for reading and posterior_sample.txt for writing
	std::fstream fin("sample.dat", std::ios::in|std::ios::binary);
	std::fstream fout("posterior_sample.txt", std::ios::out);

	fin.close();
	fout.close();
}

