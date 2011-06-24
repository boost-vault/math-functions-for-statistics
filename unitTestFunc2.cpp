// Unit test math functions for statistics proposed for TR2 - SAMPLE ONLY.

// This file is a demonstration of the proposed functions for C++ TR2
// implemented using Stephen Moshier's Cephes C library.


// Copyright Paul A. Bristow 1998-2005.
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt
// or copy at http://www.boost.org/LICENSE_1_0.txt)

// Std
#include <iostream>
	using std::cout;
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

//#include "tgmath.h" // C99 Math functions, acos, acosh ... tgamma, tfunc.
// Declares signatures for generic functions like acos(double), acos(float), acos(long double)...
// and also C names functions acosf, acosl
// (and also for complex arguments _Complex, for example cacos(), casin()...
// - but not going to attempt complex yet).

// Also useful to define template signatures for User Defined Types(UDT),
// for example, NTL by Victor Shoup,
// quad_float 106-bit precision using 128-bit
// or RR chosen-precision floating point providing 100 decimal-digit precision.
// RR and quad_float both leak memory (as advertised) but no value here so supress messages:
// BOOST_TEST_DETECT_MEMORY_LEAK 0 ?
// or post-build command line:
// "$(TargetPath)" --report_level=detailed --catch_system_errors=no --detect_memory_leak=0 --build_info=yes

// Boost
#include <boost/test/unit_test.hpp>
using boost::unit_test::test_suite;
#include <boost/test/floating_point_comparison.hpp> // Extra test tool for FP comparison.


//#ifdef _DEBUG
//  #pragma comment(lib, "i:\\lib\\libboost_unit_test_framework-vc80-mt-sgd-1_33.lib")
//#else
//  #pragma comment(lib, "i:\\boost_1_33_0\\libs\\test\\build\\msvc71_proj\\Release\\unit_test_framework.lib")
// OK 23 nOV 2005 or whne moved to I:\lib
//#endif
//#ifdef _DEBUG
//  #pragma comment(lib, "i:\\lib\\libboost_unit_test_framework-vc80-mt-sgd-1_33.lib")
//#else
//  #pragma comment(lib, "i:\\lib\\libboost_unit_test_framework-vc80-mt-s-1_33.lib")
//// But build doesn't seem to add the -s- expected.
//#endif
// But linked with these files in the property page, linker, input, additional libraries.
// and J:\Cpp\MathFuncsBoost\mathFuncLib\Debug

// C functions - 3 versions of each with suffixes for float f, and and l for long.
extern "C" double cbrt(double); // C double, global ::cbrt(9.)
extern "C" float cbrtf(float); // C float, global ::cbrt(9.F)
extern "C" double beta(double, double); //
extern "C" float betaf(float, float); //
extern "C" double acosh(double); // C99
extern "C" float acoshf(float); // C99
extern "C" double acos(double); // C99
extern "C" float acosf(float); // C99
//extern "C" double lgam(double); // log gamma Cephes - has 'wrong' name.
//extern "C" float lgamf(float); // log gamma Cephes - rename to C Std name?
// Used a macro in Cephes mconf.h when building the library to rename to C99.
extern "C" double lgamma(double); // log gamma C99 name.
extern "C" float lgammaf(float); // log gamma C99 name.
extern "C" double incbet(double, double, double); // PAB 'TR2' example.
extern "C" float incbetf(float, float, float); // PAB 'TR2' example.
extern "C" float fdtr(int, int, double); // Cephes - Fisher distribution.

#ifdef _MSC_VER
//#  pragma warning(disable: 4127) // conditional expression is constant.
//#  pragma warning(disable: 4702) // unreachable code.
//#  pragma warning(disable: 4511) // copy constructor could not be generated.
//#  pragma warning(disable: 4512) // assignment operator could not be generated.
//#  pragma warning(disable: 4521) // alignment of a member was sensitive to packing.
//#  pragma warning(disable: 4121) // alignment of a member was sensitive to packing.
//#  pragma warning(disable: 4100) // unreferenced formal parameter.
//#  pragma warning(disable: 4701) // local variable may be used without having been initialized.
//#  pragma warning(disable: 4996) // ' ' was declared deprecated.
#  pragma warning(disable: 4127) // conditional expression is constant.
#  pragma warning(disable: 4244) // conversion from 'long double' to 'double', possible loss of data
// All long double from C double produce this message at level 4 (but not 3).
#endif

// Use C function versions, for example: ::cbrt to provide C++ functions,
// specialized for float, double and long double.
// For MSVC, long double == double, so no greater accuracy possible.
// For float, using double and casting to float causes
// no loss of accuracy, but probably slower than a proper C float version.
// If float available, this is suffixed f, for example cbrtf.
// Long double cases are NOT seriously considered,
// - because MSVC long double is same as double :-(
// but Cephes DOES have long double versions, both 80-bit and 128-bit :-))

namespace std
{ // C++ Forward declarations examples for some math functions, C99/TR1.
  // To go in a mathfwd.h file?
	namespace tr1
	{ // C++ Declarations only.
		float cbrt(float);
		double cbrt(double);
		long double cbrt(long double);
	} // namespace tr1
} // namespace std

namespace std
{ // Forward declarations examples for some template functions,
	// to allow UDT to use math functions.
	namespace tr2
	{
		template <class T>	T log10(T); // Declaration to go in header.
	} // namespace tr1
} // namespace std

namespace std
{
	namespace tr1
	{ // Functions definitions added for C++ TR1 by including C99.

		double nextafter(double d, double y = numeric_limits<double>::max())
		{ // Use MS version for double.
			return _nextafter(d, y);
		} // double nextafter(double d)

		double nextbefore(double d, double y = -numeric_limits<double>::max())
		{ // Use MS version for double.
			return _nextafter(d, -y);
		} // double nextbefore(double d)

		// Cube root function.
		float cbrt(float x)
		{ // Loss of speed over C float version cbrtf?
			// But not less accurate.
			// return float(::cbrt(x));
			return ::cbrtf(x); // C99 implementation from Cephes, 32-bit IEEE 754 single precision.
		} //float cbrt(float x)
		double cbrt(double x)
		{ // Assume IEEE 754 double 64-bit precision, but warn if not.
			#if (DBL_MANT_DIG != 53)
				#pragma message ("Unknown double floating-point type!")
			#endif
			return double(::cbrt(x));
		} // double cbrt(double x)
		long double cbrt(long double x)
		{ // Loss of accuracy over C long double version cbrtl, but faster?
			#if (LDBL_MANT_DIG == 53)
						return long double(::cbrt(x)); // IEEE 754 double precision (for example MSVC).
						// Note: math.lib does not include long double cbrtl or cbrtll yet!
			#elif (LDBL_MANT_DIG == 80)
						return long double(::cbrtl(x)); // Cephes long double for IEEE 754 80-bit double extended precision.
			#elif (LDBL_MANT_DIG == 106)
						return long double(::cbrtll(x)); // Cephes 128 for IEEE 754 128-bit quadruple precision.
			#else
				#pragma message ("Unknown long double floating-point type!")
			#endif
		} // long double cbrt(long double x)

		// Log of gamma function.
		float lgamma(float x)
		{ // Log gamma function added C99.
			return float(::lgammaf(x));
		} // float lgamma(float x)
		double lgamma(double x)
		{
			return double(::lgamma(x));
		} // double lgamma(double x)
		long double lgamma(long double x)
		{
			return long double(::lgamma(x)); // lgaml
		} // long double lgamma(long double x)

		// Functions added by TR1 5.2.1 math 'special' functions.
		// [5.2.1.3] Beta function
		// beta(a, b) = gamma(a) * gamma(b) /gamma(a+b).
		// For integral values, gamma aka factorial.
		float beta(float x, float y)
		{
			return ::betaf(x, y); // global C function.
		} // float beta(float x, float y)
		double beta(double x, double y)
		{
			return ::beta(x, y);
		} // double beta(double x, double y)
		long double beta(long double x, long double y)
		{
			return (::beta(x, y));
		} // 		long double beta(long double x, long double y)
	} // namespace tr1
} // namespace std

namespace std
{
	namespace tr2
	{ // Functions 'added' by proposed 'TR2 math functions' PA Bristow WG21 N1668.
		double beta_incomplete(double a, double b, double x)
		{ // Incomplete beta function, used, in turn, by many statistical functions.
			return ::incbet(a, b, x);
		}
		float beta_incomplete(float a, float b, float x)
		{
			return ::incbetf(a, b, x);
		} // double beta_incomplete(double a, double b, double x)

		double fisher_distribution(unsigned int a, unsigned int b, double x)
		{ // Sir R A Fisher's distribution (often abbreviated to F-distribution):
			// includes the chi-squared distribution and Student's t-distribution as special cases.
			// http://en.wikipedia.org/wiki/Ronald_Fisher
			return ::fdtr(a, b, x);
		} // double fisher_distribution(unsigned int a, unsigned int b, double x)
		double log10(double x)
		{ // log to base 10.
			return ::log10(x);
		} // 		double log10(double x)

		template <class Type>
		Type log10(Type x) // Definition.
		{ // log to base 10.
			return log10(x); // Type log10
		}

	} // namespace tr2
}// namespace std

// Demonstration of a User Defined Type (UDT),
//  NTL by Victor Shoup, http://shoup.net/ntl/
#include <NTL\ZZ.h> // Arbitrary precision integer.
#include <NTL\RR.h> // Arbitrary precision floating-point.
#include <NTL\quad_float.h> // 128-bit (106-bit significand) floating point.
using NTL::quad_float;
using NTL::RR;

// Specialize ZZ, quad_float and RR for numeric_limits, as far as practicable.
// John Maddock says epsilon may be calculated from:
//   std::pow(two, 1-std::numeric_limits<T>::digits) 
using std::numeric_limits;
template <> class numeric_limits<quad_float>
 {
 public:
		static const bool is_specialized = true;
		static const int radix = 2;
		static const int digits = 106; // radix digits (bits) in significand.
		static const int digits10 = 31; // decimal digits guaranteed.
		static const int max_digits = 33; // maximum significant decimal digits.
		static const bool is_signed = true;
		static const bool is_integer = false;
		static const bool is_exact = false;
		static quad_float min() throw() { return numeric_limits<double>::min();};
		static quad_float max() throw() { return quad_float(numeric_limits<double>::max(), numeric_limits<double>::max());};
		static quad_float epsilon() throw() { return 2.4651903288156618919116517665087069677288E-32;}; // 2^(1-106) == 2^105
		static quad_float round_error() throw() { return 0.5;};
		static quad_float infinity() throw() { return quad_float(numeric_limits<double>::infinity(), numeric_limits<double>::infinity());};
		static quad_float quiet_NaN() throw() { return quad_float(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());};
		static quad_float signaling_NaN() throw() { return quad_float(numeric_limits<double>::quiet_NaN(), numeric_limits<double>::quiet_NaN());};
		static quad_float denorm_min() throw() { return quad_float(0, numeric_limits<double>::denorm_min());};
		static const int min_exponent = numeric_limits<double>::min_exponent;
		static const int min_exponent10 = numeric_limits<double>::min_exponent10;
		static const int max_exponent = numeric_limits<double>::max_exponent;
		static const int max_exponent10 = numeric_limits<double>::max_exponent10;
		static const bool has_infinity = true;
		static const bool has_quiet_NaN = true;
		static const bool has_signaling_NaN = true;
		static const float_denorm_style has_denorm = denorm_absent; // enum from <limits>
		static const bool has_denorm_loss = false;
		static const bool is_ies559 = false;
		static const bool is_bounded = true;
		static const bool is_modulo = true;
		static const bool traps = true;
		static const bool tinyness_before = true;
		static const float_round_style round_style = round_to_nearest; // enum from <limits>
 }; // class numeric_limits<quad_float>

// Need to fix precision of RR class before specialization.
// long RR::prec = 150; is in RR.cpp,
// so need a macro in RR.h to be used in RR.cpp instead of 150
#define RR_PRECISION 40

using NTL::RR;
using NTL::to_RR;
template <> class numeric_limits<RR>
 {
 public:
		static const bool is_specialized = true;
		static const int radix = 10;
		static const int digits = RR_PRECISION; // Radix digits (decimal) in significand.  ???
		static const int digits10 =RR_PRECISION -1; // Decimal digits guaranteed.
		static const int max_digits = RR_PRECISION; // Maximum significant decimal digits.
		static const bool is_signed = true;
		static const bool is_integer = false;
		static const bool is_exact = false;
		static RR min() throw() { return to_RR("0");};
		static RR max() throw() { return to_RR("99999999999999999999999999");};
		static RR epsilon() throw() { return to_RR("1e-40");};
		static RR round_error() throw() { return to_RR("0.5");};
		static RR infinity() throw() { return to_RR("9999999999999999999999999999999999999999999");};
//		static RR quiet_NaN() throw() { return RR(numeric_limits<double>::quiet_NaN())};
//		static RR signaling_NaN() throw() { return RR(numeric_limits<double>::quiet_NaN());};
//		static RR denorm_min() throw() { return RR(numeric_limits<double>::denorm_min());};
		static const int min_exponent = LONG_MIN; // numeric_limits<long>::min();
		static const int min_exponent10 = 2 + numeric_limits<long>::digits * 3010/10000; // numeric_limits<long>::digits10();
		static const int max_exponent = LONG_MAX; // numeric_limits<double>::max();
		static const int max_exponent10 = numeric_limits<double>::max_exponent10;
		static const bool has_infinity = true;
		static const bool has_quiet_NaN = true;
		static const bool has_signaling_NaN = true;
		static const float_denorm_style has_denorm = denorm_absent; // enum from <limits>
		static const bool has_denorm_loss = false;
		static const bool is_ies559 = false;
		static const bool is_bounded = true;
		static const bool is_modulo = true;
		static const bool traps = false;
		static const bool tinyness_before = false;
		static const float_round_style round_style = round_to_nearest; // enum from <limits>
 }; // class numeric_limits<RR>

using NTL::ZZ;
template <> class numeric_limits<ZZ>
{ // Defined-precision integer class.
 public:
		static const bool is_specialized = true;
		static const int radix = 2;
		static const int digits = NTL_ZZ_NBITS; // radix digits (bits).
		static const int digits10 = static_cast<int>(NTL_FRADIX); // decimal digits guaranteed.
		static const int max_digits = static_cast<int>(NTL_FRADIX); // maximum significant decimal digits.
		static const bool is_signed = true;
		static const bool is_integer = true;
		static const bool is_exact = true;
//		static ZZ min() throw() { return -(NTL_ZZ_NBITS+1);};
//		static ZZ max() throw() { return NTL_ZZ_NBITS);};
  	static const bool has_infinity = false;
		static const bool has_quiet_NaN = false;
		static const bool has_signaling_NaN = false;
		static const bool has_denorm_loss = false;
		static const bool is_ies559 = false;
		static const bool is_bounded = true;
		static const bool is_modulo = true;
		static const bool traps = false;
		static const bool tinyness_before = false;
 }; // class numeric_limits<ZZ>

union fhex
{ // Used by void outAsHex(float)
	float f;  // 32 bit real
	unsigned long l;  // 32 bit long.
};

union dhex
{ // Used by void outAsHex(double)
	double d;
	unsigned short usa[4];  // Suits Intel 8087 FP P J Plauger C Standard p 67.
};

void outAsHex(double d)
{ // displayAsHex
	dhex dd = {d};  // Initialise double dd.d = 0.0;
	// Assume Intel 8087 FP.
	// Standard C Library P J Plaguer p 67.
	char fill = cout.fill('0'); // Save.
	int fmtflags = cout.flags();
	int precision = cout.precision(2+numeric_limits<double>::digits * 3010/10000);
	cout << dec << dd.d << " == 0x" << hex << right
		<< setw(4) << dd.usa[3] // SCCC CCCC CCCC FFFF  sign & exponent.
		<< setw(4) << dd.usa[2] // CCCC
		<< setw(4) << dd.usa[1] // CCCC
		<< setw(4) << dd.usa[0] // FFFF
		;
		cout.flags(fmtflags); // Restore.
		cout.fill(fill);
		cout.precision(precision);
}  // void outAsHex(double d)

void outAsHex(quad_float l)
{ // displayAsHex
	int precision = cout.precision(2+numeric_limits<quad_float>::digits * 3010/10000); // Save.
	char fill = cout.fill('0');
	int fmtflags = cout.flags();
	cout << right << l << ' ';
	outAsHex(l.hi); cout << ' ';	outAsHex(l.lo); cout << endl;
	cout.flags(fmtflags); // Restore.
	cout.fill(fill);
	cout.precision(precision);
} //  outAsHex(quad_float l)

namespace NTL
{ // To add missing log10 function.
	quad_float log10(const quad_float& x)
	{ // log to base 10 = log(x) /log(10);
	quad_float onedivlogten;
	onedivlogten = NTL::to_quad_float("0.4342944819032518276511289189166050822944"); // Prefer * to /
	return log(x) * onedivlogten; // log(x) / log(10);
	} // quad_float log10()

	quad_float nextafter(const quad_float& x, double y = numeric_limits<double>::max())
	{ // Change by one unit in last place (least significant bit), default increases.
		// Might check that x (and y) are finite?
		// What about sign of x and y?
		double direction = (x < y) ?  +numeric_limits<double>::max() : -numeric_limits<double>::max();
		// If x < y then increment by one ulp, else decrement.
		quad_float temp = x;
		if (temp.lo == numeric_limits<double>::max())
		{ // nextafter would overflow.
			temp.lo = 0.;
			temp.hi = std::tr1::nextafter(temp.lo, direction);
		}
		else
		{ // Normal case.
			temp.lo = std::tr1::nextafter(temp.lo, direction);
		}
		return temp;
	} // quad_float nextafter(const quad_float& x)

	quad_float nextbefore(const quad_float& x)
	{ // Next lower FP value.
		quad_float temp = x;
		nextafter(temp, -numeric_limits<double>::max());
		return temp;
	} // quad_float nextbefore(const quad_float& x)
} // namespace NTL

// Implement math function test cases as free function.
// These are intended to test a few obvious, random, or exact values, a few corner cases,
// and around some points at which the algorithms are know to change,
// aiming to check for general algorithm 'sanity'
// rather than an attempt to evaluate accuracy in detail over the whole range.

// Spot values calculated using
// Stephen Moshier's Cephes DOSbox qcalc.exe 100 decimal digit accuracy,
// but displayed to 40 decimal digits (by entering "digits40")
// because this is maximum possible accuracy with foreseeable floating hardware,
// without using 'arbitrary' precision software.
// Included in Cephes 2.8 distribution http://www.moshier.net/cephes-math-28.tar.gz
// Steve Moshier's command interpreter V1.3
//
// Qlib.exe Functions are:
// h help hex acos asin atan atantwo cbrt cos cot exp expten log logten
// floor acosh asinh atanh cosh sinh tanh fac gamma lgamma jv yv ndtr
// ndtri erf erfc pdtr pdtri incbet incbetinv incgam incgaminv ellie ellik
// ellpe ellpk in gausshyp confhyp frexp ldexp polylog zeta pow sin sqrt
// tan cmp bits digits intcvts double longdouble hexinput remainder dm tm
// em take save system rem exit

void cbrt_test_function()
{ 	// Examples of cbrt, a C99 addition and hence to C++ TR1.
	BOOST_MESSAGE("Test cbrt functions.");
	BOOST_CHECK(::cbrt(-1.) == -1.); // C global
	BOOST_CHECK(::cbrt(-8.) == -2.); // Expect exact?
	BOOST_CHECK(::cbrt(0.) == 0.); // Really should be exact.
	BOOST_CHECK(::cbrt(numeric_limits<double>::max()) == 5.6438030941223622789313923870853635705819820505517385831220706618475929E102);
	BOOST_CHECK(::cbrt(27.) == 3.);
	BOOST_CHECK(::cbrt(8.) == 2.);
	BOOST_CHECK(::cbrt(2.) == 1.2599210498948731647672106072782283505702514647015079800819751121552997);
	BOOST_CHECK(std::tr1::cbrt(27.) == 3.); // C++ double version.
	BOOST_CHECK(::cbrtf(27.F) == 3.); // C global float version.
	BOOST_CHECK(std::tr1::cbrt(27.F) == 3.); // C++ float version.
	//BOOST_CHECK(std::tr1::cbrtf(27.F) == 3.); // Expected error C2039: 'cbrtf' : is not a member of 'std::tr1'
	using std::tr1::cbrt; // Convenient to avoid getting global ::cbrt by mistake.
	BOOST_CHECK(cbrt(27.) == 3.);  // C++ double version.
  BOOST_CHECK(cbrt(1.2599210498948732) == 1.08005