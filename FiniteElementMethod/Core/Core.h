#pragma once

#include "../../Matrix/Matrix.h"
#include <array>
#ifdef M_PI
#else
#define M_PI 3.14159265358979323846
#endif

namespace ComputationalMechanicsLibrary
{
	namespace FiniteElementMethod
	{
		template<typename T>
		class IElement
		{
		public:
			virtual std::vector<int>& Node() = 0;
			virtual Matrix<T>& Force() = 0;
			virtual Matrix<T>& Displacement() = 0;
			virtual Matrix<T>& Stiffness() = 0;
			virtual std::vector<int>& UnKnowForce() = 0;
			virtual std::vector<int>& UnKnowDisplacement() = 0;
		private:

		protected:

		};

		template<typename T>
		class ISolver
		{
		public:
			virtual std::vector<IElement<T>*>& Element() = 0;
			virtual Matrix<T>& Force() = 0;
			virtual Matrix<T>& Displacement() = 0;
			virtual Matrix<T>& Stiffness() = 0;
			virtual std::vector<int>& UnKnowForce() = 0;
			virtual std::vector<int>& UnKnowDisplacement() = 0;
		private:

		protected:

		};

		template<typename T>
		class ISection
		{
		public:
			virtual T Area() = 0;
			virtual T InertiaMomentY() = 0;
			virtual T InertiaMomentZ() = 0;
			virtual T PolarMoment() = 0;
		private:

		protected:

		};

		template<typename T>
		class IMaterial
		{
		public:
			virtual T ElasticityModulus() = 0;
			virtual T PoissonRatio() = 0;
			virtual T Density() = 0;
		private:

		protected:

		};

        enum DampType
        {
            Liner=1,
            FrictionCoefficient=2,
            Rayleigh=3
        };

		template<typename T>
		class IDynamicElement
		{
		public:
			virtual std::vector<int>& Node() = 0;

			virtual std::vector<Matrix<T>>& Force() = 0;

            virtual std::vector<Matrix<T>>& Displacement() = 0;
            virtual std::vector<Matrix<T>>& Velocity() = 0;
            virtual std::vector<Matrix<T>>& Acceleration() = 0;

            virtual std::vector<std::vector<int>>& UnKnowForce() = 0;
            virtual std::vector<std::vector<int>>& UnKnowDisplacement() = 0;
            virtual std::vector<std::vector<int>>& UnKnowVelocity() = 0;
            virtual std::vector<std::vector<int>>& UnKnowAcceleration() = 0;

			virtual Matrix<T>& Mass() = 0;
			virtual Matrix<T>& Damp() = 0;
			virtual Matrix<T>& Stiffness() = 0;

            virtual void DynamicReset(DampType type,int t) = 0;
		private:

		protected:

		};

		template<typename T>
		class IDynamicSolver
		{
		public:
			virtual std::vector<IDynamicElement<T>*>& Element() = 0;
			virtual std::vector<Matrix<T>>& Force() = 0;

			virtual std::vector<Matrix<T>>& Displacement() = 0;
			virtual std::vector<Matrix<T>>& Velocity() = 0;
			virtual std::vector<Matrix<T>>& Acceleration() = 0;

            virtual std::vector<std::vector<int>>& UnKnowForce() = 0;
            virtual std::vector<std::vector<int>>& UnKnowDisplacement() = 0;
            virtual std::vector<std::vector<int>>& UnKnowVelocity() = 0;
            virtual std::vector<std::vector<int>>& UnKnowAcceleration() = 0;

			virtual Matrix<T>& Mass() = 0;
			virtual Matrix<T>& Damp() = 0;
			virtual Matrix<T>& Stiffness() = 0;

			virtual T TimeInterval() = 0;
		private:

		protected:

		};

		template<typename T>
		class IArea_Section
		{
		public:
			virtual T Area() = 0;
		private:

		protected:

		};
		template<typename T>
		class IPolarInertiaMoment_Section
		{
		public:
			virtual T PolarInertiaMomentX() = 0;
			virtual T PolarInertiaMomentY() = 0;
			virtual T PolarInertiaMomentZ() = 0;
		private:

		protected:

		};
		template<typename T>
		class IInertiaMoment_Section
		{
		public:
			virtual T InertiaMomentX() = 0;
			virtual T InertiaMomentY() = 0;
			virtual T InertiaMomentZ() = 0;
		private:

		protected:

		};
		
		
		template<typename T>
		class IDynamicSection :public IArea_Section<T>,public IPolarInertiaMoment_Section<T>,public IInertiaMoment_Section<T>
		{
		public:

		private:

		protected:
		};
		
		template<typename T>
		class IDynamicMaterial:public IMaterial<T>
		{
		public:
			virtual T ShearModulus() = 0;
			virtual T ElasticityModulus() = 0;
			virtual T PoissonRatio() = 0;
			virtual T Density() = 0;
		private:

		protected:

		};


        //Function
        //
        //
        template<typename T>
        void CutZeroAddOne(ComputationalMechanicsLibrary::Matrix<T>& K,ComputationalMechanicsLibrary::Matrix<T>& D,ComputationalMechanicsLibrary::Matrix<T>& F,std::vector<int>& unKnowD,std::vector<int>& unKnowF,int length)
        {
            int newLength = 0;
            for(int f : unKnowD)
            {
                if (f == 0)
                {
                    newLength++;
                }
            }

            Matrix<T> newStiffness(K);
            Matrix<T> newForce(F);

            for (int i = 0; i<unKnowD.size(); i++)
            {
                if (unKnowD[i] != 0)
                {
                    newStiffness[i][i] = 1;
                    newForce[i][0] = D[i][0];
                    for (int j = 0; j<i; j++)
                    {
                        newStiffness[j][i] = 0;
                        newStiffness[i][j] = 0;
                    }
					for (int j = i + 1; j < length; j++)
					{
						newStiffness[j][i] = 0;
						newStiffness[i][j] = 0;
					}
				}
			}

			D = (newStiffness ^ -1) * newForce;

			F = K * D;
		}
		/** Calculate Natural Frequency with Jacobi
		 *
		 */
		template<typename T>
		Matrix<T> JacobiNaturalFrequency(Matrix<T> &KK, Matrix<T> &MM)
		{
			Matrix<T> K = KK;
			Matrix<T> M = MM;
			T a, b, c, d, α, β;
			int δ;
			Matrix<T> P(K.Row(), 1);
			Matrix<T> Q(K.Row(), 1);

			auto jacobiZero = [&](int i, int j)->void {

				if (std::abs(K[i][j]) < 1.0e-12)
				{
					return;
				}
				a = K[i][i] * M[i][j] - M[i][i] * K[i][j];
				b = K[j][j] * M[i][j] - M[j][j] * K[i][j];
				c = K[i][i] * M[j][j] - M[i][i] * K[j][j];
				if (std::abs(a) >= 1.0e-20&&std::abs(b) >= 1.0e-12)
				{
					if (c >= 0)
						δ = 1;
					else
						δ = -1;
					α =( -0.5*c + δ*std::sqrt(0.25*c*c + a*b))/a;
					β = -a*α / b;
				}
				else if (std::abs(a) < 1.0e-12&&std::abs(b) >= 1.0e-12)
				{
					α = -K[i][j] / K[j][j];
					β = 0;
				}
				else if (std::abs(a) >= 1.0e-12&&std::abs(b) < 1.0e-12)
				{
					α = 0;
					β = -K[i][j] / K[j][j];
				}
				else
				{
					α = 0;
					β = 0;
				}
				/*
				Matrix<T> E(K.Row(), K.Column());
				for (int m = 0; m < E.Row(); m++)
				{
					E[m][m] = 1;
				}
				E[i][j] = α;
				E[j][i] = β;
				*/
				T aaa = α;
				T bbb = β;

				Matrix<T> KKK = K;
				Matrix<T> MMM = M;
				//WriteMatrixToCSV("EJtemp.csv", E);
				//WriteMatrixToCSV("EJKtemp.csv", K);
				//WriteMatrixToCSV("EJTtemp.csv", E.Transport()*K);
				for (int t = 0; t < K.Row(); t++)
				{
					K[i][t] += β*K[j][t];
					K[j][t] += α*K[i][t];

					M[i][t] = MMM[i][t] + β*MMM[j][t];
					M[j][t] = MMM[j][t] + α*MMM[i][t];

				}
				for (int t = 0; t < K.Row(); t++)
				{
					K[t][i] += β*K[t][j];
					K[t][j] += α*K[t][i];

					M[t][i] += α*β*MMM[t][i] + β*MMM[t][j];
					M[t][j] += α*β*MMM[t][j] + α*MMM[t][i];
					if (std::abs(K[t][i]) <= 1.0e-12)
					{
						K[t][i] = 0;
					}
					if (std::abs(K[t][j]) <= 1.0e-12)
					{
						K[t][j] = 0;
					}
					if (std::abs(M[t][i]) <= 1.0e-12)
					{
						M[t][j] = 0;
					}
					if (std::abs(M[t][j]) <= 1.0e-12)
					{
						M[t][j] = 0;
					}
				}
				T cdefs = 0;
			};


			WriteMatrixToCSV("tempK.csv", K);
			WriteMatrixToCSV("tempM.csv", M);
			for (int k = 1; k < 1000; k++)
			{
				for (int i = 0; i < K.Column() - 1; i++)
				{
					for (int j = 0; j < i + 1; j++)
					{
						jacobiZero(j, K.Column() - 1 - i + j);
					}
				}
				WriteMatrixToCSV("jK.csv", K);
				WriteMatrixToCSV("jM.csv", M);
				for (int ii = 0; ii < K.Row(); ii++)
				{
					P[ii][0] = K[ii][ii] / M[ii][ii];
				}
				int totalEndFor = 0;
				for (int ii = 0; ii < K.Row(); ii++)
				{
					if (std::abs((P[ii][0] - Q[ii][0]) / P[ii][0]) <= 1.0e-12)
					{
						totalEndFor++;
					}
				}
				if (totalEndFor == K.Row())
				{
					break;
				}
				for (int ii = 0; ii < K.Row(); ii++)
				{
					Q[ii][0] = P[ii][0];
				}
			}

			WriteMatrixToCSV("jK.csv", K);
			WriteMatrixToCSV("jM.csv", M);
			for (int ii = 0; ii < P.Row(); ii++)
			{
				P[ii][0] = std::sqrt(std::abs(P[ii][0]));
			}
			//sort
			T iTemp;
			bool bFilish = false;
			for (int i = 0; i<P.Row() - 1; i++)
			{
				if (bFilish == true)
				{
					break;
				}
				bFilish = true;
				for (int j = P.Row() - 1; j>i; j--)
				{
					if (P[j][0]<P[j - 1][0])
					{
						iTemp = P[j - 1][0];
						P[j - 1][0] = P[j][0];
						P[j][0] = iTemp;
						bFilish = false;
					}
				}
			}

			return P;

		}

		template<typename T>
		Matrix<T> MatrixIteration(Matrix<T> &KK, Matrix<T> &MM)
		{
			Matrix<T> K = KK;
			Matrix<T> M = MM;
			
			Matrix<T> R = K^-1;
			Matrix<T> D = R*M;

			Matrix<T> A1(D.Row(), 1);
			for (int i = 0; i < A1.Row(); i++)
			{
				A1[i][0] = 1;
			}

			Matrix<T> B1 = D*A1;

			return B1;
		}

		template<typename T>
		Matrix<T> InterpolationDerivation(Matrix<T> &F)
		{
			//F.Row=2
			Matrix<T> f(F.Row(), F.Column());
			f[0][0] = 0;
			f[1][0] = F[1][0];
			for (int i = 1; i < F.Column() - 1; i++)
			{
				T leftD = (F[0][i] - F[0][i - 1]) / (F[1][i] - F[1][i - 1]);
				T rightD = (F[0][i + 1] - F[0][i]) / (F[1][i + 1] - F[1][i]);
				if(leftD*rightD<0)
				{
					f[0][i] = 0;
				}
				else
				{
					f[0][i] = leftD;
				}

				f[1][i] = F[1][i];
			}
			f[0][f.Column() - 1] = (F[0][f.Column() - 1] - F[0][f.Column() - 2]) / (F[1][f.Column() - 1] - F[1][f.Column() - 2]);

			return f;
		}

    }
}
