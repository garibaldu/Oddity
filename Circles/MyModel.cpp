#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:objects(3, 100, false, MyDistribution(-1., 1., -1., 1.))
{

}

void MyModel::fromPrior()
{
	objects.fromPrior();
}

double MyModel::perturb()
{
	double logH = 0.;

	logH += objects.perturb();

	const vector< vector<double> > components = objects.get_components();
	double rsq;
	for(size_t i=0; i<components.size(); i++)
	{
		for(size_t j=(i+1); j<components.size(); j++)
		{
			rsq = pow(components[i][1] - components[j][1], 2)
				+ pow(components[i][2] - components[j][2], 2);
			if(rsq < pow(exp(components[i][0]) + exp(components[j][0]), 2))
				return -1E250;
		}
	}		


	return logH;
}

double MyModel::logLikelihood() const
{
	double logL = 0.;

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	objects.print(out); out<<' ';
}

string MyModel::description() const
{
	return string("objects");
}

