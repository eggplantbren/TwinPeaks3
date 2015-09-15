#include "Potts.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"

#include <cmath>

using namespace DNest3;
using namespace std;

Potts::Potts()
:Model(2)
,x(50, vector<int>(50))
{

}

void Potts::from_prior()
{
	for(size_t i=0; i<x.size(); i++)
		for(size_t j=0; j<x[i].size(); j++)
			x[i][j] = randInt(num_colors);

	compute_score();
	compute_scalars();
}

void Potts::compute_score()
{
	score = 0;
	score2 = 0;

	// Coordinates of neighbouring cells
	vector<int> ii(4); vector<int> jj(4);

	for(size_t i=0; i<x.size(); i++)
	{
		for(size_t j=0; j<x[i].size(); j++)
		{
			for(int k=0; k<4; k++)
			{
				ii[k] = i;
				jj[k] = j;
			}
			// Down, up, right, left
			ii[0] = mod(i + 1, static_cast<int>(x.size()));
			ii[1] = mod(i - 1, static_cast<int>(x.size()));
			jj[2] = mod(j + 1, static_cast<int>(x[i].size()));
			jj[3] = mod(j - 1, static_cast<int>(x[i].size()));

			for(int k=0; k<4; k++)
			{
				if(x[i][j] == x[ii[k]][jj[k]])
				{
					score++;
					if(k >= 2)
						score2++;
				}
			}
		}
	}
}

#include <iostream>

double Potts::perturb()
{
	int reps = 1;
	if(rand() <= 0.5)
		reps += 1 + randInt(9);

	// Which cell is being perturbed
	int i, j;

	// Value of proposal
	int proposal;

	// Coordinates of neighbouring cells
	vector<int> ii(4); vector<int> jj(4);
	for(int z=0; z<reps; z++)
	{
		i = randInt(x.size());
		j = randInt(x[0].size());

		for(int k=0; k<4; k++)
		{
			ii[k] = i;
			jj[k] = j;
		}
		// Down, up, right, left
		ii[0] = mod(i + 1, static_cast<int>(x.size()));
		ii[1] = mod(i - 1, static_cast<int>(x.size()));
		jj[2] = mod(j + 1, static_cast<int>(x[i].size()));
		jj[3] = mod(j - 1, static_cast<int>(x[i].size()));

		// Calculate negative part of delta score
		for(int k=0; k<4; k++)
		{
			if(x[i][j] == x[ii[k]][jj[k]])
			{
				score--;
				if(k >= 2)
					score2--;
			}
		}
		// Perturb the cell
		do
		{
			proposal = randInt(num_colors);
		}while(proposal == x[i][j]);
		x[i][j] = proposal;

		// Calculate positive part of delta score
		for(int k=0; k<4; k++)
		{
			if(x[i][j] == x[ii[k]][jj[k]])
			{
				score++;
				if(k >= 2)
					score2++;
			}
		}
	}

	compute_scalars();
	return 0.;
}

void Potts::compute_scalars()
{
	scalars[0] = 0.5*score;
	scalars[1] = 0.5*score2;
}

ostream& operator << (ostream& out, const Potts& e)
{
	for(size_t i=0; i<e.x.size(); i++)
		for(size_t j=0; j<e.x[i].size(); j++)
			out<<e.x[i][j]<<' ';
	out<<e.score<<' '<<e.score2;
	return out;
}

