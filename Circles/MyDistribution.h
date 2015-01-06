#ifndef _MyDistribution_
#define _MyDistribution_

#include <Distributions/Distribution.h>

// Based on ClassicMassInf1D from RJObject
// Think of "position x" as log-period
// and mass as amplitude
class MyDistribution:public Distribution
{
	private:
		// Limits
		double x_min, x_max, y_min, y_max;
		double scale;

		// Mean and sd for logWidths
		double mu, sig;

		double perturb_parameters();

	public:
		MyDistribution(double x_min, double x_max,
				double y_min, double y_max);

		void fromPrior();

		double log_pdf(const std::vector<double>& vec) const;
		void from_uniform(std::vector<double>& vec) const;
		void to_uniform(std::vector<double>& vec) const;

		void print(std::ostream& out) const;
		static const int weight_parameter = 1;

};

#endif

