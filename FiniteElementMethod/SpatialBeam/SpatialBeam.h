#pragma once
#include "../Core/Core.h"
#define M_PI 3.14159265358979323846

namespace ComputationalMechanicsLibrary
{
	namespace FiniteElementMethod
	{
		namespace SpatialBeam
		{
			template<typename T>
			class SpatialBeamElement :public IDynamicElement<T>
			{
			public:
				SpatialBeamElement(int _i,int _j, IDynamicSection<T>& _section, double _length, double _angle, IDynamicMaterial<T>& _Material, Matrix<T>& _force, Matrix<T>& _displacement, Matrix<T> _velocity,Matrix<T> _acceleration, std::vector<int>& _unKnowForce, std::vector<int>& _unKnowDisplacement,std::vector<int>& _unKnowVelocity,std::vector<int>& _unKnowAcceleration)
				{
					this->node = { _i,_j };
					this->sectionArea = _section.Area();
					this->polarInertiaMomentX = _section.PolarInertiaMomentX();
					this->polarInertiaMomentY = _section.PolarInertiaMomentY();
					this->polarInertiaMomentZ = _section.PolarInertiaMomentZ();
					this->inertiaMomentX = _section.InertiaMomentX();
					this->inertiaMomentY = _section.InertiaMomentY();
					this->inertiaMomentZ = _section.InertiaMomentZ();
					this->length = _length;
					this->angle = _angle;
					this->shearModulus = _Material.ShearModulus();
					this->elasticityModulus = _Material.ElasticityModulus();
					this->poissonRatio = _Material.PoissonRatio();
					this->density = _Material.Density();
					this->force = _force;
					this->displacement = _displacement;
					this->velocity = _velocity;
					this->acceleration = _acceleration;
					this->unKnowForce = _unKnowForce;
					this->unKnowDisplacement = _unKnowDisplacement;
					this->unKnowVelocity = _unKnowVelocity;
					this->unKnowAcceleration = _unKnowAcceleration;

					this->CalculateMass();
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
				Matrix<T>& Velocity()
				{
					return this->velocity;
				}
				Matrix<T>& Acceleration()
				{
					return this->acceleration;
				}

				Matrix<T>& Mass()
				{
					return this->mass;
				}
				Matrix<T>& Damp()
				{
					return this->damp;
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
				std::vector<int>& UnKnowVelocity()
				{
					return this->unKnowVelocity;
				}
				std::vector<int>& UnKnowAcceleration()
				{
					return this->unKnowAcceleration;
				}
			private:
				std::vector<int> node;

				Matrix<T> force;

				Matrix<T> displacement;
				Matrix<T> velocity;
				Matrix<T> acceleration;

				Matrix<T> mass;
				Matrix<T> damp;
				Matrix<T> stiffness;

				std::vector<int>& unKnowForce;
				std::vector<int>& unKnowDisplacement;
				std::vector<int>& unKnowVelocity;
				std::vector<int>& unKnowAcceleration;

				T shearModulus;
				T elasticityModulus;
				T poissonRatio;
				T density;
				T sectionArea;
				T polarInertiaMomentX;
				T polarInertiaMomentY;
				T polarInertiaMomentZ;
				T inertiaMomentX;
				T inertiaMomentY;
				T inertiaMomentZ;
				T length;
				T angle;

				void CalculateMass()
				{
					T l = this->length;
					T IxA = this->inertiaMomentX / this->sectionArea;
					this->mass = [[140, 0, 0, 0, 0, 0, 70, 0, 0, 0, 0, 0],
						[0, 156, 0, 0, 0, 22 * l, 0, 54, 0, 0, 0, -13 * l],
						[0, 0, 156, 0, -22 * l, 0, 0, 0, 54, 0, 13 * l, 0],
						[0, 0, 0, 40 * IxA, 0, 0, 0, 0, 0, 70 * IxA, 0, 0],
						[0, 0, -22 * l, 0, 4 * l*l, 0, 0, 0, -13 * l, 0, -3 * l*l, 0],
						[0, 22 * l, 0, 0, 0, 4 * l*l, 0, 13 * l, 0, 0, 0, -3 * l*l],
						[70, 0, 0, 0, 0, 0, 140, 0, 0, 0, 0, 0],
						[0, 54, 0, 0, 0, 13 * l, 0, 156, 0, 0, 0, -22 * l],
						[0, 0, 54, 0, -13 * l, 0, 0, 0, 156, 0, 22 * l, 0],
						[0, 0, 0, 70 * IxA, 0, 0, 0, 0, 0, 140 * IxA, 0, 0],
						[0, 0, 13 * l, 0, -3 * l*l, 0, 0, 0, 22 * l, 0, 4 * l*l, 0],
						[0, -13 * l, 0, 0, 0, -3 * l*l, 0, -22 * l, 0, 0, 0, 4 * l*l]];
					this->mass = (this->density*this->sectionArea*this->length / 420)* this->mass;
				}
				void CalculateStiffness()
				{
					T& dix = this->displacement[0];
					T& djx = this->displacement[6];
					T& ¦Èix = this->displacement[3];
					T& ¦Èjx = this->displacement[9];
					T& ¦Í = this->poissonRatio;
					T& G = this->shearModulus;
					T& E = this->elasticityModulus;
					T& A = this->sectionArea;
					T& l = this->length;
					T& Ix = this->inertiaMomentX;
					T& Iy = this->inertiaMomentY;
					T& Iz = this->inertiaMomentZ;
					T& Jx = this->polarInertiaMomentX;
					T& Jy = this->polarInertiaMomentY;
					T& Jz = this->polarInertiaMomentZ;

					Matrix<T> Kl = [[E*A / l, 0, 0, 0, 0, 0, -E*A / l, 0, 0, 0, 0, 0],
						[0, 12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l, 0, -12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l],
						[0, 0, 12 * E*Iy / l / l / l, 1, -6 * E*Iy / l / l, 0, 0, 0, -12 * E*Iy / l / l / l, 0, -6 * E*Iy / l / l, 0],
						[0, 0, 0, G*Jx / l, 0, 0, 0, 0, 0, -G*Jx / l, 0, 0],
						[0, 0, -6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0],
						[0, 6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l, 0, 0, -6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l],
						[-E*A / l, 0, 0, 0, 0, 0, E*A / l, 0, 0, 0, 0, 0],
						[0, -12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l, 0, 12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l],
						[0, 0, -12 * E*Iy / l / l / l, 0, 6 * E*Iz / l / l, 0, 0, 12 * E*Iy / l / l / l, 0, 6 * E*Iy / l / l, 0],
						[0, 0, 0, -G*Jx / l, 0, 0, 0, 0, 0, G*Jx / l, 0, 0],
						[0, 0, -6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0],
						[0, 6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l, 0, -6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l]];
					
					Matrix<T> KNB = [[3 / 2, 0, 0, 0, 0, 0, -3 / 2, 0, 0, 0, 0, 0],
						[0, 6 / 5, 0, 0, 0, l / 10, 0, -6 / 5, 0, 0, 0, l / 10],
						[0, 0, 6 / 5, 0, -l / 10, 0, 0, 0, -6 / 5, 0, -l / 10, 0],
						[0, 0, 0, Ix / A, 0, 0, 0, 0, 0, -Ix / A, 0, 0],
						[0, 0, -l / 10, 0, 2 * l*l / 15, 0, 0, 0, l / 10, 0, -l*l / 30, 0],
						[l / 10, 0, 0, 0, 0, 2 * l*l / 15, -l / 10, 0, 0, 0, 0, -l*l / 30],
						[-3 / 2, 0, 0, 0, 0, 0, 3 / 2, 0, 0, 0, 0, 0],
						[0, -6 / 5, 0, 0, 0, -l / 10, 0, 6 / 5, 0, 0, 0, -l / 10],
						[0, 0, -6 / 5, 0, l / 10, 0, 0, 0, 6 / 5, 0, l / 10, 0],
						[0, 0, 0, -Ix / A, 0, 0, 0, 0, 0, ix / A, 0, 0],
						[0, 0, -l / 10, 0, -l*l / 30, 0, 0, 0, l / 10, 0, 2 * l*l / 15, 0],
						[l / 10, 0, 0, 0, 0, -l*l / 30, -l / 10, 0, 0, 0, 0, 2 * l*l / 15]];
					KNB = E*A*(djx - dix) / l / l*KNB;

					Matrix<T> KNT = [[0, 0, 0, ¦Í / 2, 0, 0, 0, 0, 0, -¦Í / 2, 0, 0],
						[0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0],
						[0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1],
						[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						[0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -l / 2],
						[0, 0, -1, 0, 0, 0, 0, 0, 1, 0, l / 2, 0],
						[0, 0, 0, -¦Í / 2, 0, 0, 0, 0, 0, ¦Í / 2, 0, 0],
						[0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0],
						[0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1],
						[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
						[0, 1, 0, 0, 0, l / 2, 0, -1, 0, 0, 0, 0],
						[0, 0, 1, 0, -¦Í / 2, 0, 0, 0, -1, 0, 0, 0]];
					KNT = (1 + 2 * ¦Í)*G*Ix*(¦Èjx - ¦Èix) / l / l * KNT;

					this->stiffness = Kl + KNB + KNT;
				}
			protected:

			};

			template<typename T>
			class SpatialBeamSolver :public IDynamicSolver<T>
			{
			public:
				SpatialBeamSolver(std::vector<IDynamicElement<T>*>& _element)
				{
					//assemb integral matrix
					this->element = _element;
					int length = this->Length();

					this->mass=Matrix<T>(length, length);
					this->damp = Matrix<T>(length, length);
					this->stiffness = Matrix<T>(length, length);

					this->force.push_back(Matrix<T>(length, 1));
					this->displacement.push_back(Matrix<T>(length, 1));
					this->velocity.push_back(Matrix<T>(length, 1));
					this->acceleration.push_back(Matrix<T>(length, 1));

					this->unKnowForce.resize(length);
					this->unKnowDisplacement.resize(length);
					this->unKnowVelocity.resize(length);
					this->unKnowAcceleration.resize(length);

					for (IDynamicElement<T>* e : this->element)
					{
						int ei = e->Node()[0];
						int ej = e->Node()[1];

						int ni = ei, nj = ei;
						for (int i = 0; i < e->Stiffness().Row(); i++)
						{
							if (i >= 3)
							{
								ni = ej;
							}
							else
							{
								ni = ei;
							}

							this->force[0][(ni - 1) * 3 + i % 3][0] += e->Force()[i][0];
							this->unKnowForce[(ni - 1) * 3 + i % 3] = e->UnKnowForce()[i];
							this->displacement[0][(ni - 1) * 3 + i % 3][0] += e->Displacement()[i][0];
							this->unKnowDisplacement[(ni - 1) * 3 + i % 3] = e->UnKnowDisplacement()[i];
							this->velocity[0][(ni - 1) * 3 + i % 3][0] += e->Velocity()[i][0];
							this->unKnowVelocity[(ni - 1) * 3 + i % 3] = e->UnKnowVelocity()[i];
							this->acceleration[0][(ni - 1) * 3 + i % 3][0] += e->Acceleration()[i][0];
							this->unKnowAcceleration[(ni - 1) * 3 + i % 3] = e->UnKnowAcceleration()[i];


							for (int j = 0; j < e->Stiffness().Column(); j++)
							{
								if (j >= 3)
								{
									nj = ej;
								}
								else
								{
									nj = ei;
								}

								this->mass[(ni - 1) * 3 + i % 3][(nj - 1) * 3 + j % 3] += e->Mass()[i][j];
								this->damp[(ni - 1) * 3 + i % 3][(nj - 1) * 3 + j % 3] += e->Damp()[i][j];
								this->stiffness[(ni - 1) * 3 + i % 3][(nj - 1) * 3 + j % 3] += e->Stiffness()[i][j];
							}
						}


					}

				}

				std::vector<IDynamicElement<T>*>& Element()
				{
					return this->element;
				}
				std::vector<Matrix<T>>& Force()
				{
					return this->force;
				}

				std::vector<Matrix<T>>& Displacement()
				{
					return this->displacement;
				}
				std::vector<Matrix<T>>& Velocity()
				{
					return this->velocity;
				}
				std::vector<Matrix<T>>& Acceleration()
				{
					return this->acceleration;
				}

				Matrix<T>& Mass()
				{
					return this->mass;
				}
				Matrix<T>& Damp()
				{
					return this->damp;
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
				std::vector<int>& UnKnowVelocity()
				{
					return this->unKnowVelocity;
				}
				std::vector<int>& UnKnowAcceleration()
				{
					return this->unKnowAcceleration;
				}

				int Length()
				{
					//consecutive node number : 1,2,3,4,...
					int all = 0;
					for (IDynamicElement<T>* e : this->element)
					{
						for (int i : e->Node())
						{
							if (i > all)
							{
								all = i;
							}
						}
					}

					return all * 6;
				}
			private:
				std::vector<IElement<T>*> element;
				std::vector<Matrix<T>> force;

				std::vector<Matrix<T>> displacement;
				std::vector<Matrix<T>> velocity;
				std::vector<Matrix<T>> acceleration;

				Matrix<T> mass;
				Matrix<T> damp;
				Matrix<T> stiffness;

				std::vector<int> unKnowForce;
				std::vector<int> unKnowDisplacement;
				std::vector<int> unKnowVelocity;
				std::vector<int> unKnowAcceleration;

			protected:

			};

			template<typename T>
			class SpatialBeamMaterial :public IDynamicMaterial<T>
			{
			public:
				SpatialBeamMaterial(T _elasticityModulus, T _poissonRatio, T _density)
				{
					this->elasticityModulus = _elasticityModulus;
					this->poissonRatio = _poissonRatio;
					this->density = _density;
				}

				T ShearModulus()
				{
					return this->elasticityModulus / (2 * (1 + this->poissonRatio));
				}
				T ElasticityModulus()
				{
					return this->elasticityModulus;
				}
				T PoissonRatio()
				{
					return this->poissonRatio;
				}
				T Density()
				{
					return this->density;
				}
			private:
				T elasticityModulus;
				T poissonRatio;
				T density;
			protected:

			};

			template<typename T>
			class SpatialBeamSection :public IDynamicSection<T>
			{
			public:
				SpatialBeamSection(T _diameter = 0)
				{
					this->Diameter = _diameter;
				}

				T Area()
				{
					return M_PI* std::pow(this->Radius(), 2);
				}
				
				T PolarInertiaMomentX()
				{
					return M_PI * std::pow(this->Radius(), 4) / 2;
				}
				T PolarInertiaMomentY()
				{
					return M_PI * std::pow(this->Radius(), 4) / 2;
				}
				T PolarInertiaMomentZ()
				{
					return M_PI * std::pow(this->Radius(), 4) / 2;
				}
				T InertiaMomentX()
				{
					return M_PI * std::pow(this->Radius(), 4) / 2;
				}
				T InertiaMomentY()
				{
					return M_PI * std::pow(this->Radius(), 4) / 4;
				}
				T InertiaMomentZ()
				{
					return M_PI * std::pow(this->Radius(), 4) / 4;
				}

				T Radius()
				{
					return this->Diameter / 2;
				}

				T Diameter;
			private:

			protected:

			};
		}
	}
}