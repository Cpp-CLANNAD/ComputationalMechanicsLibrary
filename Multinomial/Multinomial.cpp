#include "Multinomial.h"

int ComputationalMechanicsLibrary::Factorial(int x)
{
	if (x <= 1)
	{
		return 1;
	}
	else
	{
		return Factorial(x - 1)*x;
	}
}

int ComputationalMechanicsLibrary::Combinatorial(int n, int m)
{
	return Factorial(m) / (Factorial(m-n)*Factorial(n));
}
