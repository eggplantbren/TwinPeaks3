#include <iostream>
#include <valarray>

using namespace std;

extern "C"
{
	void __junk_MOD_increment(double* x, int* M, int* N);
}


int main()
{
	valarray<double> x{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

	int M{2}, N{3};
	__junk_MOD_increment(&x[0], &M, &N);

	for(const double& xx: x)
		cout<<xx<<' ';
	cout<<endl;

	return 0;
}

