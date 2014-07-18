#include <iostream>
#include <ctime>
#include <cstdlib>
#include <RandomNumberGenerator.h>
#include "Sampler.h"
#include "Models/Potts.h"

using namespace std;
using namespace DNest3;

int main()
{
	char choice;
	do
	{
		cout<<"# [O]verwrite or [A]ppend output files? ";
		cin>>choice;
		choice = tolower(choice);
	}while(choice != 'o' && choice != 'a');
	if(choice == 'o')
	{
		int result = system("rm logw_thinned.txt logw.txt sample.txt scalars_thinned.txt scalars.txt");
		if(result == 0)
			cout<<"# Files removed."<<endl;
		else
			cout<<"# Files didn't exist."<<endl;
	}

	// Initialise RNG and seed with time
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	Sampler<Potts> s(100, 500, 2000);
	s.initialise();

	for(int i=0; i<200000; i++)
		s.update();

	return 0;
}

