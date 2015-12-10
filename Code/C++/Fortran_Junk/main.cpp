#include <iostream>
#include <valarray>

using namespace std;

extern "C"
{
	void __junk_MOD_increment(double* x, int* N);
}

void increment(valarray<double>& x)
{
	int N = x.size();
	__junk_MOD_increment(&x[0], &N);
}

int main()
{
	valarray<double> x{1.0, 2.0, 3.0, 4.0, 5.0};
	increment(x);

	for(const double& xx: x)
		cout<<xx<<' ';
	cout<<endl;

	return 0;
}

