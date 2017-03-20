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
                    for (int j = i + 1; j<length; j++)
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
			for (int k = 1; k < 1000; k++)
			{
				for (int i = 0; i < K.Row(); i++)
				{
					for (int j = i + 1; j < K.Row(); j++)
					{
						a = K[i][i] * M[i][j] - M[i][i] * K[i][j];
						b = K[j][j] * M[i][j] - M[j][j] * K[i][j];
						c = K[i][i] * M[j][j] - M[i][i] * K[j][j];
						if (c >= 0)
							δ = 1;
						else
							δ = -1;
						d = 0.5*c + δ*sqrt((0.5*c)*(0.5*c) + a*b);
						α = d <= 1e-20 ? 0 : b / d;
						β = d <= 1e-20 ? 0 : -a / d;
						Matrix<T> E(K.Row(), K.Column());
						for (int m = 0; m < E.Row(); m++)
						{
							E[m][m] = 1;
						}
						E[i][j] = α;
						E[j][i] = β;

						//WriteMatrixToCSV("EJtemp.csv", E);
						//WriteMatrixToCSV("EJKtemp.csv", K);
						//WriteMatrixToCSV("EJTtemp.csv", E.Transport()*K);
						K = E.Transport()*K*E;
						M = E.Transport()*M*E;
						for (int m = 0; m < K.Row(); m++)
						{
							for (int n = 0; n < K.Column(); n++)
							{
								if (std::abs(K[m][n]) <= 1e-20)
								{
									K[m][n] = 0;
								}
								if (std::abs(M[m][n]) <= 1e-20)
								{
									M[m][n] = 0;
								}
							}
						}
					}
				}
				//WriteMatrixToCSV("tempK.csv", K);
				//WriteMatrixToCSV("tempM.csv", M);
				for (int ii = 0; ii < K.Row(); ii++)
				{
					P[ii][0] = K[ii][ii] / M[ii][ii];
				}
				int totalEndFor = 0;
				for (int ii = 0; ii < K.Row(); ii++)
				{
					if (std::abs((P[ii][0] - Q[ii][0]) / P[ii][0]) <= 1e-20)
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

			//WriteMatrixToCSV("aK.csv", K);
			//WriteMatrixToCSV("aM.csv", M);
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

    }
}
