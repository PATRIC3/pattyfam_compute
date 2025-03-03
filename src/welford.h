#ifndef _welford_h
#define _welford_h

#include <cmath>

/*
 * Welford's algorithm for online computation of mean & variance.
 * https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Online_algorithm
 */

class Welford
{
public:
    int count_ = 0;
    double mean_ = 0.0;
    double M2_ = 0.0;

    void update(double new_value) {
	count_++;
	double delta = new_value - mean_;
	mean_ += delta / (double) count_;
	double delta2 = new_value - mean_;
	M2_ += delta * delta2;
    }

    void finalize(double &mean, double &variance, double &sample_variance) {
	if (count_ < 2)
	{
	    mean = variance = sample_variance = std::nan("");
	}
	else
	{
	    mean = mean_;
	    variance = M2_ / static_cast<double>(count_);
	    sample_variance = M2_ / static_cast<double>(count_ - 1);
	}
    }
};


#endif
