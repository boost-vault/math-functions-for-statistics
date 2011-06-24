// math2.cpp
// A few samples of proposed  math functions for TR2
// Implemented using Stephen Moshier's Cephes C functions library.

extern "C" double incbet(double, double, double); // PAB 'TR2' example.
extern "C" float incbetf(float, float, float); // PAB 'TR2' example.
extern "C" float fdtr(int, int, double); // Cephes - Fisher distribution.
extern "C" float fdtrc(int, int, double); // Cephes - Fisher distribution complement.


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
			return ::fdtrc(a, b, x);
		} // double fisher_distribution(unsigned int a, unsigned int b, double x)
	} // namespace tr2
}// namespace std
