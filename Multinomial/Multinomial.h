#pragma once
#include <cmath>
#include <functional>
#include "../../ComputationalMechanicsLibrary/Matrix/Matrix.h"

namespace ComputationalMechanicsLibrary
{
	template<typename T>
	class Multinomial
	{
	public:
		T* coefficients;
		int n;
		Multinomial(int _n=1)
		{
			this->coefficients = nullptr;
			this->n = _n;

			this->coefficients = new T[_n];
			for (int i = 0; i < _n; i++)
			{
				this->coefficients[i] = 0;
			}
		}
		Multinomial(Multinomial &m)
		{
			this->coefficients = nullptr;
			this->n = m.n;

			this->coefficients = new T[m.n];
			for (int i = 0; i < m.n; i++)
			{
				this->coefficients[i] = m.coefficients[i];
			}
		}

		T Solve(T x)
		{
			T result = 0;
			for (int i = 0; i < this->n; i++)
			{
				result += this->coefficients[i] * std::pow(x, i);
			}

			return result;
		}

		Multinomial<T> Differentiate()
		{
			Multinomial<T> result(this->n == 1 ? 1 : this->n - 1);

			for (int i = 0; i < result.n; i++)
			{
				result.coefficients[i] = this->coefficients[i + 1] * (i + 1);
			}

			return result;
		}

		Multinomial<T> Multiply(Multinomial<T> &b)
		{
			Multinomial<T> result(this->n + b.n - 1);
			for (int i = 0; i < this->n; i++)
			{
				for (int j = 0; j < b.n; j++)
				{
					result.coefficients[i + j] += this->coefficients[i] * b.coefficients[j];
				}
			}

			return result;
		}

		Multinomial<T> Add(Multinomial<T> &b)
		{
			Multinomial<T> result((b.n > this->n )? b.n : this->n);

			for (int i = 0; i < result.n; i++)
			{
				if (i < b.n)
				{
					result.coefficients[i] += b.coefficients[i];
				}
				if (i < this->n)
				{
					result.coefficients[i] += this->coefficients[i];
				}
			}

			return result;
		}

		Multinomial<T> X0Expand(T x0)
		{
			Multinomial<T> result(1);
			Multinomial<T> *axs = new Multinomial<T>[this->n];
			for (int i = 0; i<this->n; i++)
			{
				axs[i] = Multinomial<T>(i + 1);
				for (int j = 0; j < i + 1; j++)
				{
					(axs[i]).coefficients[j] += this->coefficients[i]*std::pow(x0, i - j)*Combinatorial(j, i);
				}
				result = result.Add(axs[i]);
			}
			delete[] axs;
			return result;
		}

		Matrix<T> ToMatrix()
		{
			Matrix<T> result(1, this->n);
			for (int i = 0; i < this->n; i++)
			{
				result[0][i] = this->coefficients[i];
			}

			return result;
		}

		~Multinomial()
		{
			//delete[] this->coefficients;
		}
	};

	template<typename T>
	Multinomial<T> Taylor(std::function<T(T)> func, T x0 = 0, int N = 3)
	{
		Multinomial<T> result(N);
		T deltaX = 0.0001;
		//differentiate
		std::function<T(T, int, std::function<T(T)>, T)> dfn;
		dfn = [&dfn](T x, int n, std::function<T(T)> f, T dx)->T {
			if (n == 0)
			{
				return f(x);
			}
			else
			{
				return (dfn(x + dx, n - 1, f, dx) - dfn(x, n - 1, f, dx)) / dx;
			}
		};

		//factorial
		std::function<T(T)> nnn;
		nnn = [&nnn](T n)->T {
			if (n <= 1)
			{
				return 1;
			}
			else
			{
				return nnn(n - 1)*n;
			}
		};

		for (int i = 0; i < N; i++)
		{
			result.coefficients[i] = dfn(x0, i, func, deltaX) / nnn(i);
		}

		return result;
	}
	template<typename T>
	Multinomial<T> Taylor(std::function<T(T)> funcX, std::function<T(T)> funcY, T t0 = 0, int N = 3)
	{
		Multinomial<T> result(N);
		T deltaT = 0.0001;

		for (int i = 0; i < N; i++)
		{
			result.coefficients[i] = NumericalDifferentiation(t0, i, funcY, deltaT) / NumericalDifferentiation(t0, i, funcX, deltaT) / Factorial(i);
			deltaT = 0.0001;
		}

		return result;
	}

	int Factorial(int x);

	int Combinatorial(int n, int m);

	template<typename T>
	T NumericalDifferentiation(T x, int n, std::function < T(T)> f, T dx)
	{
		if (n == 0)
		{
			return f(x);
		}
		else
		{
			return (NumericalDifferentiation(x + dx, n - 1, f, dx) - NumericalDifferentiation(x, n - 1, f, dx)) / dx;
		}
	}

	template<typename T>
	std::function<T(T)> Hermite(Matrix<T> xy)
	{
		std::function<T(T)> fx = [&](T x) {
			return x;
		};


		return fx;
	}

	template<typename T>
	class ScatterFunction
	{
	public:
		Matrix<T> xy;

		ScatterFunction(Matrix<T> _xy)
		{
			this->xy = _xy;

			BuildFunction();
		}

		void BuildFunction()
		{
			T df0 = (this->xy[1][2] - this->xy[1][0]) / (this->xy[0][2] - this->xy[0][0]);
			T dfn = (this->xy[1][this->xy.Column() - 1] - this->xy[1][this->xy.Column() - 3]) / (this->xy[0][this->xy.Column() - 1] - this->xy[0][this->xy.Column() - 3]);

			Matrix<T> h(1, this->xy.Column());
			for (int i = 1; i < h.Column(); i++)
			{
				h[0][i] = this->xy[0][i] - this->xy[0][i - 1];
			}
			Matrix<T> miu(1, this->xy.Column());
			Matrix<T> lambda(1, this->xy.Column());
			for (int i = 1; i < this->xy.Column() - 1; i++)
			{
				miu[0][i] = h[0][i] / (h[0][i + 1] + h[0][i]);
				lambda[0][i] = 1 - miu[0][i];
			}


			Matrix<T> A(this->xy.Column() - 1, this->xy.Column() - 1);
			Matrix<T> B(this->xy.Column() - 1, 1);
			A[0][0] = 2;
			A[0][1] = 1;
			B[0][0] = 6 * (this->difference(this->xy[0][0], this->xy[0][1]) - df0) / (h[0][1]);

			for (int i = 1; i < A.Row() - 1; i++)
			{
				A[i][i - 1] = miu[0][i];
				A[i][i] = 2;
				A[i][i + 1] = lambda[0][i];
				B[i][0] = 6 * this->difference(this->xy[0][i - 1], this->xy[0][i], this->xy[0][i + 1]);
			}

			A[A.Row() - 1][A.Column() - 2] = 1;
			A[A.Row() - 1][A.Column() - 1] = 2;
			B[B.Row() - 1][0] = 6 * (dfn - this->difference(this->xy[0][this->xy.Column() - 2], this->xy[0][this->xy.Column() - 1])) / h[0][h.Column() - 1];

			this->M = (A^-1) * B;
		}

		T Solve(T x)
		{
			for (int i = 0; i < this->xy.Column() - 1; i++)
			{
				if (x >= this->xy[0][i] && x <= this->xy[0][i + 1])
				{
					T hi = this->xy[0][i + 1] - this->xy[0][i];
					return this->M[i][0] * std::pow(this->xy[0][i + 1] - x, 3) / (6 * hi) +
						this->M[i + 1][0] * std::pow(x - this->xy[0][i], 3) / (6 * hi) +
						(this->xy[1][i] / hi - M[i][0] * hi / 6)*(this->xy[0][i + 1] - x) +
						(this->xy[1][i + 1] / hi - M[i + 1][0] * hi / 6)*(x - this->xy[0][i]);
				}
			}
		}
	private:
		Matrix<T> M;
		//²î·Ö
		T difference(T x0,T x1)
		{
			return (differenceSolve(x1) - differenceSolve(x0)) / (x1 - x0);
		}
		T difference(T x0,T x1,T x2)
		{
			return (difference(x2, x1) - difference(x1, x0)) / (x2 - x0);
		}

		T differenceSolve(T x)
		{
			for (int i = 0; i < this->xy.Column(); i++)
			{
				if (this->xy[0][i] == x)
				{
					return this->xy[1][i];
				}
			}
		}
	};

	template<typename T>
	Multinomial<T> LeastSquaresMethodToMultinomial(Matrix<T> xy,int k)
	{
		std::function<T(int)> sumxPow = [&](int n)
		{
			T sxn = 0;
			for (int i = 0; i < xy.Column(); i++)
			{
				sxn += std::pow(xy[0][i], n);
			}
			return sxn;
		};
		std::function<T(int)> sumxPowy = [&](int n)
		{
			T sxny = 0;
			for (int i = 0; i < xy.Column(); i++)
			{
				sxny += std::pow(xy[0][i], n)*xy[1][i];
			}
			return sxny;
		};
		Multinomial<T> result(k + 1);

		Matrix<T> A(k + 1, k + 1);
		Matrix<T> B(k + 1, 1);
		for (int i = 0; i < A.Row(); i++)
		{
			for (int j = 0; j < A.Column(); j++)
			{
				A[i][j] = sumxPow(i + j);
			}

			B[i][0] = sumxPowy(i);
		}
		WriteMatrixToCSV("A.csv", A);
		WriteMatrixToCSV("B.csv", B);
		Matrix<T> a = (A^-1)*B;
		for (int i = 0; i < k + 1; i++)
		{
			result.coefficients[i] = a[i][0];
		}

		return result;
	}

	template<typename T>
	Matrix<T> Differentiation(Matrix<T> &f)
	{
		Matrix<T> result = f;

		result[1][0] = (f[1][1] - f[1][0]) / (f[0][1] - f[0][0]);

		for (int i = 1; i < result.Column() - 1; i++)
		{
			result[1][i] = ((f[1][i + 1] - f[1][i]) / (f[0][i + 1] - f[0][i]) + (f[1][i] - f[1][i - 1]) / (f[0][i] - f[0][i - 1])) / 2;
		}

		result[1][result.Column() - 1] = (f[1][f.Column() - 1] - f[1][f.Column() - 2]) / (f[0][f.Column() - 1] - f[0][f.Column() - 2]);

		return result;
	}
}
