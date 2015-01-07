#include "MyDistribution.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include <cmath>
#include <gsl/gsl_cdf.h>

using namespace DNest3;

MyDistribution::MyDistribution(double x_min, double x_max,
				double y_min, double y_max)
:x_min(x_min)
,x_max(x_max)
,y_min(y_min)
,y_max(y_max)
,scale(sqrt((x_max-x_min)*(y_max-y_min)))
{

}

void MyDistribution::fromPrior()
{
	mu = log(1E-3*scale) + log(1E3)*randomU();
	sig = 0.1 + 2.9*randomU();
}

double MyDistribution::perturb_parameters()
{
	double logH = 0.;

	int which = randInt(2);

	if(which == 0)
	{
		mu += log(1E3)*randh();
		wrap(mu, log(1E-3*scale), log(scale));
	}
	if(which == 1)
	{
		sig += 2.9*randh();
		wrap(sig, 0.1, 3.);
	}

	return logH;
}

// vec[0] = logWidth
// vec[1] = xc
// vec[2] = yc

double MyDistribution::log_pdf(const std::vector<double>& vec) const
{
	double logP = 0.;

	logP += -log(sig) - 0.5*pow((vec[0] - mu)/sig, 2);
	if(vec[1] < x_min || vec[1] > x_max
	|| vec[2] < y_min || vec[2] > y_max)
		logP = -1E300;

	return logP;
}

void MyDistribution::from_uniform(std::vector<double>& vec) const
{
	vec[0] = mu + sig*gsl_cdf_ugaussian_Pinv(vec[0]);
	vec[1] = x_min + (x_max - x_min)*vec[1];
	vec[2] = y_min + (y_max - y_min)*vec[2];
}

void MyDistribution::to_uniform(std::vector<double>& vec) const
{
	vec[0] = gsl_cdf_ugaussian_P((vec[0] - mu)/sig);
	vec[1] = (vec[1] - x_min)/(x_max - x_min);
	vec[2] = (vec[2] - y_min)/(y_max - y_min);
}

void MyDistribution::print(std::ostream& out) const
{
	out<<mu<<' '<<sig<<' ';
}

