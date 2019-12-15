#include <iostream>
#include <random>
#include "random.h"

using namespace std;


NormalDistribution::NormalDistribution(double mean, double _std) {
	time_t seed  = time(nullptr);
	engine       = new default_random_engine(seed);
	distribution = new normal_distribution<double>(mean, _std);
}


double NormalDistribution::read() {
	return (*distribution)(*engine);
}