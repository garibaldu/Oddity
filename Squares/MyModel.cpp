#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>

using namespace std;
using namespace DNest3;

MyModel::MyModel()
:objects(3, 1000, false, MyDistribution(-1., 1., -1., 1.))
{

}

void MyModel::fromPrior()
{
	objects.fromPrior();
	calculate_overlap();
}

double MyModel::perturb()
{
	double logH = 0.;

	int old_overlap = overlap;

	logH += objects.perturb();
	calculate_overlap();

	if(overlap > old_overlap)
		logH = -1E200;

	return logH;
}

void MyModel::calculate_overlap()
{
	overlap = 0;

	const vector< vector<double> > components = objects.get_components();
	vector<double> widths(components.size());
	for(size_t i=0; i<components.size(); i++)
		widths[i] = exp(components[i][0]);

	double diffx, diffy;
	for(size_t i=0; i<components.size(); i++)
	{
		for(size_t j=(i+1); j<components.size(); j++)
		{
			diffx = fabs(components[i][1] - components[j][1]);
			diffy = fabs(components[i][2] - components[j][2]);
			if(diffx < (widths[i] + widths[j]) && (diffy < (widths[i] + widths[j])))
				overlap++;
		}
	}
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

