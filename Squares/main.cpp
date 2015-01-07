#include <iostream>
#include "Start.h"
#include "MyModel.h"

using namespace std;
using namespace DNest3;

int main(int argc, char** argv)
{
	MTSampler<MyModel> sampler = setup_mt<MyModel>(argc, argv);
	sampler.run();

	return 0;
}

