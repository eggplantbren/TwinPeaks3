#include <vector>
#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include "RandomNumberGenerator.h"

using namespace std;
using namespace DNest3;

class Data
{
	private:
		vector<double> logw, run_id;
		vector< vector<double> > scalars;
		int N;

	public:
		Data()
		{ }

		void load()
		{
			logw.clear(); run_id.clear(); scalars.clear();

			double temp;
			fstream fin("logw.txt", ios::in);
			int k = 0;
			while(fin >> temp)
			{
				logw.push_back(temp);
				if(k == 0)
					run_id.push_back(0);
				else if(logw.back() == -1.)
					run_id.push_back(run_id.back() + 1);
				else
					run_id.push_back(run_id.back());
				k++;
			}
			fin.close();

			fin.open("scalars.txt", ios::in);
			double temp2;
			while(fin >> temp && fin >> temp2)
			{
				vector<double> vec(2);
				vec[0] = temp; vec[1] = temp2;
				scalars.push_back(vec);
			}
			fin.close();

			if(scalars.size() < logw.size())
			{
				logw.erase(logw.begin() + scalars.size(),
						logw.end());
				run_id.erase(run_id.begin() + scalars.size(),
						run_id.end());
			}
			if(logw.size() < scalars.size())
			{
				scalars.erase(scalars.begin() + logw.size(),
						scalars.end());
			}

			assert(logw.size() == scalars.size());
			assert(logw.size() == run_id.size());
			N = logw.size();
			cout<<"# Loaded "<<N<<" points."<<endl;
	}

	friend class Assignment;
};


/************************************************************************/

class Assignment
{
	private:
		const Data& data;
		vector<double> logX;

	public:
		Assignment(const Data& data)
		:data(data)
		{
		}

		void initialise()
		{
			logX.assign(data.N, 0.);

			for(int i=0; i<data.N; i++)
			{
				if(data.logw[i] == -1.)
					logX[i] = log(randomU());
				else
					logX[i] = logX[i-1] + log(randomU());			
			}
		}

		int update_all()
		{
			int total = 0;
			for(int i=0; i<data.N; i++)
			{
				total += update(i);
				cout<<"."<<flush;
			}
			cout<<endl;
			return total;
		}

		// Gibbs sample one value
		int update(int i)
		{
			double proposal;

			// Range of values for proposal
			double lower = -1E300;
			double upper = 0.;

			// Lower limit -- next point in same run (if it exists)
  			if( (i != (data.N-1)) &&
					(data.run_id[i+1] == data.run_id[i]))
    				lower = logX[i+1];
			// Upper limit -- previous point in same run (if it exists)
  			if( (i != 0) &&
					(data.run_id[i-1] == data.run_id[i]))
    				upper = logX[i-1];

			// If lower limit exists, generate uniformly between limits
			if(lower != -1E300)
				proposal = lower + (upper - lower)*randomU();
			else // otherwise use exponential distribution
				proposal = upper + log(randomU());

			// Measure inconsistency
			int inconsistency_old = 0;
			int inconsistency_new = 0;

			bool outside, inside;
			for(int j=0; j<data.N; j++)
			{
				if(data.run_id[j] != data.run_id[i])
				{
					// See if distribution j is outside,
					// inside, or unknown
					outside = true;
					inside  = true;
					for(int k=0; k<2; k++)
					{
						if(data.scalars[j][k] > data.scalars[i][k])
							outside = false;
						if(data.scalars[j][k] < data.scalars[i][k])
							inside = false;
					}

					if(outside && (logX[j] < logX[i]))
						inconsistency_old++;
					if(inside  && (logX[j] > logX[i]))
						inconsistency_old++;

					if(outside && (logX[j] < proposal))
						inconsistency_new++;
					if(inside  && (logX[j] > proposal))
						inconsistency_new++;
				}
			}

			int inconsistency = inconsistency_old;
			if(inconsistency_new <= inconsistency_old)
			{
				logX[i] = proposal;
				inconsistency = inconsistency_new;
			}
			return inconsistency;
		}

		void save()
		{
			fstream fout("logX.txt", ios::out);
			for(int i=0; i<data.N; i++)
				fout<<logX[i]<<endl;
			fout.close();
		}

};


int main()
{
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	Data data;
	data.load();

	Assignment assignment(data);
	assignment.initialise();

	int inconsistency;
	int when_finished = -1;
	int i=0;
	while(true)
	{
		inconsistency = assignment.update_all();
		assignment.save();
		if(when_finished == -1 && inconsistency == 0)
			when_finished = i;
		i++;
		cout<<i<<' '<<inconsistency<<endl;
		if(i > (1 + 1.2*when_finished))
			break;
	}
	return 0;
}

