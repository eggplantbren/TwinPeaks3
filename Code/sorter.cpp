#include <iostream>
#include <vector>
#include <fstream>
#include <RandomNumberGenerator.h>

using namespace DNest3;
using namespace std;

/************************************************************/

class Point
{
	private:
		std::vector<double> scalars;
		int run, step_within_run;

	public:
		Point(double S1, double S2, int run, int step_within_run)
		:scalars(2)
		,run(run)
		,step_within_run(step_within_run)
		{
			scalars[0] = S1;
			scalars[1] = S2;
		}

		friend bool operator < (const Point& p1, const Point& p2);
		friend ostream& operator << (ostream& out, const Point& p);
};

ostream& operator << (ostream& out, const Point& p)
{
	for(size_t i=0; i<p.scalars.size(); i++)
		out<<p.scalars[i]<<' ';
	out<<p.run<<' '<<p.step_within_run;
	return out;
}

bool operator < (const Point& p1, const Point& p2)
{
	if(p1.run == p2.run)
		return p1.step_within_run < p2.step_within_run;

	for(size_t i=0; i<p1.scalars.size(); i++)
		if(p1.scalars[i] > p2.scalars[i])
			return false;
	return true;
}

/************************************************************/

using namespace std;

int main()
{
	// Initialise RNG and seed with time
	RandomNumberGenerator::initialise_instance();
	RandomNumberGenerator::get_instance().set_seed(time(0));

	// Size of the data we're dealing with
	int steps = 200;

	// Load all the points
	fstream fin("scalars.txt", ios::in);
	vector<Point> points;
	double S1, S2;
	int k = 0;
	while(fin>>S1 && fin>>S2)
	{
		points.push_back(Point(S1, S2, k/steps, k%steps));
		k++;
	}
	fin.close();
	cout<<"# Loaded "<<points.size()<<" points. ";

	// Number of runs
	int runs = (int)(points.size())/steps;

	cout<<"That's "<<runs<<" runs."<<endl;

	sort(points.begin(), points.end());

	fstream fout("sorted.txt", ios::out);
	for(size_t i=0; i<points.size(); i++)
		fout<<points[i]<<endl;
	fout<<endl;

	return 0;
}

