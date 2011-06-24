namespace std
{
	namespace tr2
	{ // Functions 'added' by proposed 'TR2 math functions' PA Bristow WG21 N1668.
		double beta_incomplete(double a, double b, double x);
		 // Incomplete beta function, used, in turn, by many statistical functions.
		float beta_incomplete(float a, float b, float x);

		double fisher_distribution(unsigned int a, unsigned int b, double x);
		// Sir R A Fisher's distribution (often abbreviated to F-distribution):
		// includes the chi-squared distribution and Student's t-distribution as special cases.
		// http://en.wikipedia.org/wiki/Ronald_Fisher
		// double fisher_distribution(unsigned int a, unsigned int b, double x);
	} // namespace tr2
}// namespace std
