#pragma once

#include "../../Matrix/Matrix.h"
#include <array>

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

		template<typename T>
		class IDynamicElement :public IElement<T>
		{
		public:
			virtual std::vector<int>& Node() = 0;

			virtual Matrix<T>& Force() = 0;

			virtual Matrix<T>& Displacement() = 0;
			virtual Matrix<T>& Velocity() = 0;
			virtual Matrix<T>& Acceleration() = 0;

			virtual Matrix<T>& Mass() = 0;
			virtual Matrix<T>& Damp() = 0;
			virtual Matrix<T>& Stiffness() = 0;

			virtual std::vector<int>& UnKnowForce() = 0;
			virtual std::vector<int>& UnKnowDisplacement() = 0;
			virtual std::vector<int>& UnKnowVelocity() = 0;
			virtual std::vector<int>& UnKnowAcceleration() = 0;
		private:

		protected:

		};

		template<typename T>
		class IDynamicSolver:public ISolver<T>
		{
		public:
			virtual std::vector<IDynamicElement<T>*>& Element() = 0;
			virtual std::vector<Matrix<T>>& Force() = 0;

			virtual std::vector<Matrix<T>>& Displacement() = 0;
			virtual std::vector<Matrix<T>>& Velocity() = 0;
			virtual std::vector<Matrix<T>>& Acceleration() = 0;

			virtual Matrix<T>& Mass() = 0;
			virtual Matrix<T>& Damp() = 0;
			virtual Matrix<T>& Stiffness() = 0;

			virtual std::vector<int>& UnKnowForce() = 0;
			virtual std::vector<int>& UnKnowDisplacement() = 0;
			virtual std::vector<int>& UnKnowVelocity() = 0;
			virtual std::vector<int>& UnKnowAcceleration() = 0;
		private:

		protected:

		};

		template<typename T>
		class IArea_Section :public ISection<T>
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
		class IDynamicSection :public IArea_Section<T>, IPolarInertiaMoment_Section<T>, IInertiaMoment_Section<T>
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
	}
}