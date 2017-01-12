#pragma once
#include "../Core/Core.h"
#define M_PI 3.14159265358979323846

namespace ComputationalMechanicsLibrary
{
	namespace FiniteElementMethod
	{
		namespace PlaneRigidFrame
		{
			/// <summary>
			/// 2 node plane rigid frame element
			/// </summary>
			template<typename T>
			class PlaneRigidFrameElement :public IElement<T>
			{
			public:

				PlaneRigidFrameElement(int _i, int _j, ISection<T>& _section, double _length, double _angle, IMaterial<T>& _Material, Matrix<T>& _force, Matrix<T>& _displacement, std::vector<int>& _unKnowForce, std::vector<int>& _unKnowDisplacement)
				{
					this->node = { _i,_j };
					this->sectionArea = _section.Area();
					this->inertiaMoment = _section.InertiaMomentZ();
					this->length = _length;
					this->angle = _angle;
					this->elasticityModulus = _Material.ElasticityModulus();
					this->force = _force;
					this->displacement = _displacement;
					this->unKnowForce = _unKnowForce;
					this->unKnowDisplacement = _unKnowDisplacement;

					this->CalculateStiffness();
				}

				std::vector<int>& Node()
				{
					return this->node;
				}
				Matrix<T>& Force()
				{
					return this->force;
				}
				Matrix<T>& Displacement()
				{
					return this->displacement;
				}
				Matrix<T>& Stiffness()
				{
					return this->stiffness;
				}
				std::vector<int>& UnKnowForce()
				{
					return this->unKnowForce;
				}
				std::vector<int>& UnKnowDisplacement()
				{
					return this->unKnowDisplacement;
				}
			private:
				std::vector<int> node;
				Matrix<T> force;
				Matrix<T> displacement;
				Matrix<T> stiffness;
				std::vector<int> unKnowForce;
				std::vector<int> unKnowDisplacement;

				T elasticityModulus;
				T sectionArea;
				T inertiaMoment;
				T length;
				T angle;

				void CalculateStiffness()
				{
					T cc = std::cos(this->angle);
					T ss = std::sin(this->angle);
					Matrix<T> Transport = { { cc, ss, 0, 0, 0, 0 },{ -ss, cc, 0, 0, 0, 0 },{ 0, 0, 1, 0, 0, 0 },{ 0, 0, 0, cc, ss, 0 },{ 0, 0, 0, -ss, cc, 0 },{ 0, 0, 0, 0, 0, 1 } };

					T s = this->elasticityModulus * this->inertiaMoment / (this->length * this->length * this->length);
					T a = this->elasticityModulus * this->sectionArea / this->length;
					T b = s * 12;
					T c = s * 6 * this->length;
					T d = s * 4 * this->length * this->length;
					this->stiffness = { { a, 0, 0, -a, 0, 0 },{ 0, b, c, 0, -b, c },{ 0, c, d, 0, -c, d / 2 },{ -a, 0, 0, a, 0, 0 },{ 0, -b, -c, 0, b, -c },{ 0, c, d / 2, 0, -c, d } };
					this->stiffness = (Transport.Transport()) * (this->stiffness) * Transport;
				}
			protected:

			};
		
			template<typename T>
			class CircularSection :public ISection<T>
			{
			public:
				CircularSection(T _diameter=0)
				{
					this->Diameter = _diameter;
				}

				T Area()
				{
					return M_PI* std::pow(this->Radius(), 2);
				}
				T InertiaMomentY()
				{
					return M_PI* std::pow(this->Radius(), 4) / 4;
				}
				T InertiaMomentZ()
				{
					return M_PI* std::pow(this->Radius(), 4) / 4;
				}
				T PolarMoment()
				{
					return M_PI * std::pow(this->Radius(), 4) / 2;
				}

				T Radius()
				{
					return this->Diameter / 2;
				}

				T Diameter;
			private:

			protected:

			};

			template<typename T>
			class PlaneRigidFrameMaterial:public IMaterial<T>
			{
			public:
				PlaneRigidFrameMaterial(T _elasticityModulus, T _poissonRatio)
				{
					this->elasticityModulus = _elasticityModulus;
					this->poissonRatio = _poissonRatio;
				}

				T ElasticityModulus()
				{
					return this->elasticityModulus;
				}
				T PoissonRatio()
				{
					return this->poissonRatio;
				}
			private:
				T elasticityModulus;
				T poissonRatio;
			protected:

			};
		}
	}
}

