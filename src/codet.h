#ifndef _codet_h
#define _codet_h

/**
 * Helper to compute coefficient of determination.
 */

class CoefficentOfDetermination
{
public:
    CoefficentOfDetermination(double len)
	: len_(len)
	, sum_(0.0)
	, sum_res_(0.0) {
    }
    void insert(double val) {
	val /= len_;
	vals_.push_back(val);
	sum_ += val;
//	double ei = val - len_;
	double ei = val - 1.0;
	sum_res_ += ei * ei;
	sum_dev_ += abs(ei);
    }
    double deviation() {
	return sum_dev_ / static_cast<double>(vals_.size()) / len_;
    }
    double compute() {
	double sum_tot = 0.0;
	double m = sum_ / static_cast<double>(vals_.size());
	for (auto v: vals_)
	{
	    sum_tot += (v - m) * (v - m);
	}
	double R2 = 1.0 - sum_res_ / sum_tot;
	// std::cerr << "m=" << m << " sum_res=" << sum_res_ << " sum_tot=" << sum_tot << " R2=" << R2 << "\n";
	return R2;
    }

private:
    double len_;
    double sum_;
    double sum_res_;
    double sum_dev_;
    std::vector<double> vals_;
};

#endif




