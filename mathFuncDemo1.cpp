// Demonstration of math functions for statistics proposed for TR2.

// This file is a simple demonstration of one of the proposed
// mathematical functions (incomplete beta) for C++ TR2
// implemented using Stephen Moshier's Cephes C library.

// Copyright Paul A. Bristow 1998-2005.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Std
#include <iostream>
	using std::cout;
	using std::cerr;
	using std::cin;
	using std::endl;
	using std::dec;
	using std::hex;
	using std::right;
	using std::boolalpha;
#include <iomanip>
	using std::setprecision;
	using std::setw;
#include <string>
	using std::string;
#include <limits>
	using std::numeric_limits;
#include <cmath>
	double sqrt(double);

// Sample math functions provided by the proposed std::tr2 library.
// Declared in math2.cpp, defined in math2.hpp, wrapping C library functions.

// double beta_incomplete(double a, double b, double x);
// Incomplete beta function.
// http://en.wikipedia.org/wiki/Incomplete_beta_function

// double fisher_distribution(unsigned int a, unsigned int b, double x);
// Sir R A Fisher's distribution (often abbreviated to F-distribution):
// includes the chi-squared distribution and Student's t-distribution as special cases.
// http://en.wikipedia.org/wiki/Ronald_Fisher

#include "math2.hpp"
  using std::tr2::beta_incomplete;
  using std::tr2::fisher_distribution;

#ifdef _MSC_VER
#  pragma warning(disable: 4127) // conditional expression is constant.
#endif

// J. C. Miller and J. N. Miller, Statistics For Analytical Chemistry,
// 3rd edition, Ellis Horwoord ISBN 0 13030990-7
// Section 3.3 comparison of means of two samples
// pages 56 -57.
// Determination of tin in foodstuffs refluxing with acid for 30 and 75 minutes.
// tin found (mg/kg)
// Anaytical Methods Committee, analyst (1983) vol 109, page 109.
// Does the mean amount of tin found differ significantly for the two refluxing times?

double tin1 [] = {55., 57., 59., 56., 56., 59.};
double tin2 [] = {57., 55., 58., 59., 59., 59.};

template <typename Type>
inline Type sqr(const Type& a)
{ // Notationally convenient for a^2 (but overflow possible!)
	return (Type)a * a;
}

void meanvariance(double data[], size_t n, double* mean, double* variance)
{ // Compute mean and variance of first n items in data array.
	*mean = 0.;
	for (size_t j = 0; j < n; j++)
	{
		*mean += data[j];
	}
	*mean /= n; // Final mean.
	double s = 0.; // Sum of differences from mean.
	*variance = 0.; // Sum of squared differences.
	for (size_t j = 0; j < n; j++)
	{
		double d = data[j]-(*mean); // Difference from mean.
		s += d;
		*variance += d * d;
	}
	*variance = (*variance - s * s/n)/(n-1);
} // void meanvariance(double data[], size_t n, double *mean, double *variance)

int main()
{
	cout << "Demo of comparing means using incomplete beta math function.";
# if defined(__FILE__) && defined(__TIMESTAMP__)
	cout << "  " << __FILE__ << ' ' << __TIMESTAMP__ << ' '<< _MSC_FULL_VER << endl;
# endif

	const size_t n1 = sizeof(tin1)/sizeof(double);
	const size_t n2 = sizeof(tin2)/sizeof(double);
	if (n1 != n2)
	{
		cout << "number of determinations differ,"
			"method 1  " << n1 << ", method 2 " << n2 << endl;
	}
	else
	{
		cout << n1 << " determinations of lead concentration by both methods." << endl;
	}
	cout.precision(3); // To suit precision of data items (about 1%).

	double mean_1;
	double variance_1;
	double mean_2;
	double variance_2;

	meanvariance(tin1, n1, &mean_1, &variance_1);
	cout << "Mean of 1 method is " << mean_1 << endl;
	cout << "Variance of 1 method is " << variance_1 << endl;
	meanvariance(tin2, n2, &mean_2, &variance_2);
	cout << "Mean of 2 method is " << mean_2 << endl;
	cout << "Variance of 2 method is " << variance_2 << endl;
	cout << "Difference of means is " << mean_1 - mean_2 << endl;
	cout << "Difference of variance  is " << variance_1 - variance_2 << endl;
	cout.precision(2);  // To suit precision of statistical parameters with only 10 degrees of freedom.
	double t = (mean_1 - mean_2)/sqrt(variance_1/n1 + variance_2/n2);
	cout << "Student's t is " << t << endl;
	double df = sqr(variance_1/n1 + variance_2/n2) /(sqr(variance_1/n1) / (n1-1) + sqr(variance_2/n2)/ (n2-1));
	cout << "Degrees of freedom is " << df << endl;

	using std::tr2::beta_incomplete; // Needed to avoid ambiguity with C99 function definition.
	// error C2668: 'std::tr2::beta_incomplete' : ambiguous call to overloaded function
	// could be 'double std::tr2::beta_incomplete(double,double,double)'
	// or 'double beta_incomplete(double,double,double)'
	// while trying to match the argument list '(double, double, double)'
	double prob = beta_incomplete(0.5 * df, 0.5, df/(df + t * t));
	cout << "Probability of 1 method being same as 2 method is "  << prob << endl;
	if (prob > 0.05)
	{
		cout << "We can be >95% confident that longer boiling does not change the lead released." << endl;
	}

	double f; // Variance ratio.
	size_t df1;
	size_t df2;
	if (variance_1 > variance_2)
	{ // ensure F is > 1.
		f = variance_1/variance_2;
		df1 = n1-1;
		df2 = n2-1;
	}
	else
	{
		f = variance_2/variance_1;
		df1 = n2-1;
		df2 = n1-1;
	}
	cout << "Fisher F is " << f << endl;
	prob = 2. * fisher_distribution(n1-1, n2-1, f); // Small probability means variances are different.
	// if (prob >1) prob = 1. - prob;
	cout << "Probability that methods are equally precise is " << prob  << endl;
	if (prob > 0.05)
	{
		cout << "We can be >95% confident that the methods are equally precise." << endl;
	}
	cout << endl;
	return 0;
}  // int main()


/*

Output is
------ Build started: Project: MathFuncDemo1, Configuration: Release Win32 ------
Compiling...
mathFuncDemo1.cpp
MSVC++ compiler Version 8.0
Linking...
Generating code
Finished generating code
Autorun "j:\Cpp\MathFuncsBoost\release\MathFuncDemo1.exe"
Demo of comparing means using incomplete beta math function.
.\mathFuncDemo1.cpp Mon Nov 28 19:04:37 2005 140050727


6 determinations of lead concentration by both methods.
Mean of 1 method is 57
Variance of 1 method is 2.8
Mean of 2 method is 57.8
Variance of 2 method is 2.57
Difference of means is -0.833
Difference of variance  is 0.233
Student's t is -0.88
Degrees of freedom is 10
Probability of 1 method being same as 2 method is 0.4
We can be >95% confident that longer boiling does not change the lead released.
Fisher F is 1.1
Probability that methods are equally precise is 0.93
We can be >95% confident that the methods are equally precise.


Build log was saved at "file://j:\Cpp\MathFuncsBoost\MathFuncDemo1\Release\BuildLog.htm"
MathFuncDemo1 - 0 error(s), 0 warning(s)
========== Build: 1 succeeded, 0 failed, 0 up-to-date, 0 skipped ==========
*/