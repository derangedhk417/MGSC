#include <iostream>
#include <random>
#include <chrono>

using namespace std;

class NormalDistribution {
public:
	default_random_engine       * engine;
	normal_distribution<double> * distribution;
	NormalDistribution(double mean, double _std);
	double read();
};