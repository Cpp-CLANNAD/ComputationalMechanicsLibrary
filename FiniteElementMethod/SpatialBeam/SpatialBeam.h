#pragma once
#include "../Core/Core.h"
#include "../../Multinomial/Multinomial.h"

#include <functional>

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
				SpatialBeamElement(int _i, int _j, IDynamicSection<T>& _section, double _length,double _azimuthAngle, double _wellAngle, IDynamicMaterial<T>& _Material,
					std::vector<Matrix<T>>& _force, std::vector<Matrix<T>>& _displacement, std::vector<Matrix<T>>& _velocity, std::vector<Matrix<T>>& _acceleration,
					std::vector<std::vector<int>>& _unKnowForce, std::vector<std::vector<int>>& _unKnowDisplacement, std::vector<std::vector<int>>& _unKnowVelocity, std::vector<std::vector<int>>& _unKnowAcceleration,
					DampType _type = DampType::Rayleigh)
				{
					this->type = _type;
					this->node = { _i,_j };
					this->sectionArea = _section.Area();
					this->polarInertiaMomentX = _section.PolarInertiaMomentX();
					this->polarInertiaMomentY = _section.PolarInertiaMomentY();
					this->polarInertiaMomentZ = _section.PolarInertiaMomentZ();
					this->inertiaMomentX = _section.InertiaMomentX();
					this->inertiaMomentY = _section.InertiaMomentY();
					this->inertiaMomentZ = _section.InertiaMomentZ();
					this->length = _length;
					this->wellAngle = _wellAngle;
					this->azimuthAngle = _azimuthAngle;
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

					this->DynamicResetFunctions;

					this->CalculateTransformation();
					this->CalculateMass();
					this->CalculateStiffness();
					this->CalculateDamp(this->type);

					for (Matrix<T> &f : this->force)
					{
						f = this->transformation.Transport()*f;
					}

					this->stiffness = this->transformation.Transport()*this->stiffness*this->transformation;

				}

				std::vector<int>& Node()
				{
					return this->node;
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

				std::vector<std::vector<int>>& UnKnowForce()
				{
					return this->unKnowForce;
				}
				std::vector<std::vector<int>>& UnKnowDisplacement()
				{
					return this->unKnowDisplacement;
				}
				std::vector<std::vector<int>>& UnKnowVelocity()
				{
					return this->unKnowVelocity;
				}
				std::vector<std::vector<int>>& UnKnowAcceleration()
				{
					return this->unKnowAcceleration;
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

				T Length()
				{
					return this->length;
				}
				T Density()
				{
					return this->density;
				}
				T SectionArea()
				{
					return this->sectionArea;
				}
				T WellAngle()
				{
					return this->wellAngle;
				}
				T AzimuthAngle()
				{
					return this->azimuthAngle;
				}


				std::vector<std::function<void(SpatialBeamElement<T>*,int t)>> DynamicResetFunctions;

				void DynamicReset(DampType type, int t)
				{
					/*
					this->force[t] = this->transformation*this->force[t];
					this->displacement[t] = this->transformation*this->displacement[t];
					this->velocity[t] = this->transformation*this->velocity[t];
					this->acceleration[t] = this->transformation*this->acceleration[t];*/
					
					for (Matrix<T> &f : this->force)
					{
						f = this->transformation*f;
					}
					for (Matrix<T> &d : this->displacement)
					{
						d = this->transformation*d;
					}
					for (Matrix<T> &v : this->velocity)
					{
						v = this->transformation*v;
					}
					for (Matrix<T> &a : this->acceleration)
					{
						a = this->transformation*a;
					}
					
					this->CalculateStiffness();
					this->CalculateDamp(type);

					for (auto DynamicResetFunction : this->DynamicResetFunctions)
					{
						DynamicResetFunction(this, t);
					}
					
					/*
					this->force[t] = this->transformation.Transport()*this->force[t];
					this->displacement[t] = this->transformation.Transport()*this->displacement[t];
					this->velocity[t] = this->transformation.Transport()*this->velocity[t];
					this->acceleration[t] = this->transformation.Transport()*this->acceleration[t];*/
					
					for (Matrix<T> &f : this->force)
					{
						f = this->transformation.Transport()*f;
					}
					for (Matrix<T> &d : this->displacement)
					{
						d = this->transformation.Transport()*d;
					}
					for (Matrix<T> &v : this->velocity)
					{
						v = this->transformation.Transport()*v;
					}
					for (Matrix<T> &a : this->acceleration)
					{
						a = this->transformation.Transport()*a;
					}
					
				}
			private:
				DampType type;

				std::vector<int> node;

				std::vector<Matrix<T>> force;

				std::vector<Matrix<T>> displacement;
				std::vector<Matrix<T>> velocity;
				std::vector<Matrix<T>> acceleration;

				std::vector<std::vector<int>> unKnowForce;
				std::vector<std::vector<int>> unKnowDisplacement;
				std::vector<std::vector<int>> unKnowVelocity;
				std::vector<std::vector<int>> unKnowAcceleration;

				Matrix<T> mass;
				Matrix<T> damp;
				Matrix<T> stiffness;

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
				T wellAngle;
				T azimuthAngle;

				Matrix<T> transformation;

				void CalculateMass()
				{
					T l = this->length;
					T IxA = this->inertiaMomentX / this->sectionArea;
					this->mass = { { 140, 0, 0, 0, 0, 0, 70, 0, 0, 0, 0, 0 },
					{ 0, 156, 0, 0, 0, 22 * l, 0, 54, 0, 0, 0, -13 * l },
					{ 0, 0, 156, 0, -22 * l, 0, 0, 0, 54, 0, 13 * l, 0 },
					{ 0, 0, 0, 140 * IxA, 0, 0, 0, 0, 0, 70 * IxA, 0, 0 },
					{ 0, 0, -22 * l, 0, 4 * l*l, 0, 0, 0, -13 * l, 0, -3 * l*l, 0 },
					{ 0, 22 * l, 0, 0, 0, 4 * l*l, 0, 13 * l, 0, 0, 0, -3 * l*l },
					{ 70, 0, 0, 0, 0, 0, 140, 0, 0, 0, 0, 0 },
					{ 0, 54, 0, 0, 0, 13 * l, 0, 156, 0, 0, 0, -22 * l },
					{ 0, 0, 54, 0, -13 * l, 0, 0, 0, 156, 0, 22 * l, 0 },
					{ 0, 0, 0, 70 * IxA, 0, 0, 0, 0, 0, 140 * IxA, 0, 0 },
					{ 0, 0, 13 * l, 0, -3 * l*l, 0, 0, 0, 22 * l, 0, 4 * l*l, 0 },
					{ 0, -13 * l, 0, 0, 0, -3 * l*l, 0, -22 * l, 0, 0, 0, 4 * l*l } };
					this->mass = (this->density*this->sectionArea*this->length / 420)* this->mass;
				}
				void CalculateStiffness()
				{
					T dix = this->displacement[0][0][0];
					T djx = this->displacement[0][6][0];
					T θix = this->displacement[0][3][0];
					T θjx = this->displacement[0][9][0];
					T ν = this->poissonRatio;
					T G = this->shearModulus;
					T E = this->elasticityModulus;
					T A = this->sectionArea;
					T l = this->length;
					T Ix = this->inertiaMomentX;
					T Iy = this->inertiaMomentY;
					T Iz = this->inertiaMomentZ;
					T Jx = this->polarInertiaMomentX;
					T Jy = this->polarInertiaMomentY;
					T Jz = this->polarInertiaMomentZ;

					Matrix<T> Kl = { { E*A / l, 0, 0, 0, 0, 0, -E*A / l, 0, 0, 0, 0, 0 },
					{ 0, 12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l, 0, -12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l },
					{ 0, 0, 12 * E*Iy / l / l / l,0, -6 * E*Iy / l / l, 0, 0, 0, -12 * E*Iy / l / l / l, 0, -6 * E*Iy / l / l, 0 },
					{ 0, 0, 0, G*Jx / l, 0, 0, 0, 0, 0, -G*Jx / l, 0, 0 },
					{ 0, 0, -6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0 },
					{ 0, 6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l, 0, 0, -6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l },
					{ -E*A / l, 0, 0, 0, 0, 0, E*A / l, 0, 0, 0, 0, 0 },
					{ 0, -12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l, 0, 12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l },
					{ 0, 0, -12 * E*Iy / l / l / l, 0, 6 * E*Iz / l / l, 0, 0, 12 * E*Iy / l / l / l, 0, 6 * E*Iy / l / l, 0 },
					{ 0, 0, 0, -G*Jx / l, 0, 0, 0, 0, 0, G*Jx / l, 0, 0 },
					{ 0, 0, -6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0 },
					{ 0, 6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l, 0, -6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l } };

					Matrix<T> KNB = { { 3.0 / 2, 0, 0, 0, 0, 0, -3.0 / 2, 0, 0, 0, 0, 0 },
					{ 0, 6.0 / 5, 0, 0, 0, l / 10, 0, -6.0 / 5, 0, 0, 0, l / 10 },
					{ 0, 0, 6.0 / 5, 0, -l / 10, 0, 0, 0, -6.0 / 5, 0, -l / 10, 0 },
					{ 0, 0, 0, Ix / A, 0, 0, 0, 0, 0, -Ix / A, 0, 0 },
					{ 0, 0, -l / 10, 0, 2 * l*l / 15, 0, 0, 0, l / 10, 0, -l*l / 30, 0 },
					{ l / 10, 0, 0, 0, 0, 2 * l*l / 15, -l / 10, 0, 0, 0, 0, -l*l / 30 },
					{ -3.0 / 2, 0, 0, 0, 0, 0, 3.0 / 2, 0, 0, 0, 0, 0 },
					{ 0, -6.0 / 5, 0, 0, 0, -l / 10, 0, 6.0 / 5, 0, 0, 0, -l / 10 },
					{ 0, 0, -6.0 / 5, 0, l / 10, 0, 0, 0, 6.0 / 5, 0, l / 10, 0 },
					{ 0, 0, 0, -Ix / A, 0, 0, 0, 0, 0, Ix / A, 0, 0 },
					{ 0, 0, -l / 10, 0, -l*l / 30, 0, 0, 0, l / 10, 0, 2 * l*l / 15, 0 },
					{ l / 10, 0, 0, 0, 0, -l*l / 30, -l / 10, 0, 0, 0, 0, 2 * l*l / 15 } };
					KNB = E*A*(djx - dix) / l / l*KNB;

					Matrix<T> KNT = { { 0, 0, 0, ν / 2, 0, 0, 0, 0, 0, -ν / 2, 0, 0 },
					{ 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 1.0, 0 },
					{ 0, 0, 0, 0, 0, -1.0, 0, 0, 0, 0, 0, 1.0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, -1.0, 0, 0, 0, 0, 0, 1.0, 0, 0, 0, -l / 2 },
					{ 0, 0, -1.0, 0, 0, 0, 0, 0, 1.0, 0, l / 2, 0 },
					{ 0, 0, 0, -ν / 2, 0, 0, 0, 0, 0, ν / 2, 0, 0 },
					{ 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, -1.0, 0 },
					{ 0, 0, 0, 0, 0, 1.0, 0, 0, 0, 0, 0, -1.0 },
					{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
					{ 0, 1.0, 0, 0, 0, l / 2, 0, -1.0, 0, 0, 0, 0 },
					{ 0, 0, 1.0, 0, -ν / 2, 0, 0, 0, -1.0, 0, 0, 0 } };
					KNT = (1 + 2 * ν)*G*Ix*(θjx - θix) / l / l * KNT;

					this->stiffness = Kl + KNB + KNT;
				}

				void CalculateDamp(DampType type = DampType::Rayleigh)
				{
					switch (type) {
					case DampType::Liner:
						this->LinerDamp();
						break;
					case DampType::FrictionCoefficient:
						this->FrictionCoefficientDamp();
						break;
					case DampType::Rayleigh:
						this->RayleighDamp();
						break;
					default:
						this->RayleighDamp();
						break;
					}
				}

				void LinerDamp()
				{
					T l = this->length;
					this->damp = { { (l / 3),0,0,0,0,0,(l / 6),0,0,0,0,0 },
					{ 0,-(16.0 / 15) + 2.0 / l + (4 * l) / 5,0,0,0,1.0 / 10 + l / 30 + l*l / 30,0,8.0 / 15 - 6.0 / (5 * l) - (3 * l) / 10,0,0,0,1 - (2 * l) / 3 + l*l / 60 },
					{ 0,0,6.0 / (5 * l) + (13 * l) / 35,0,-(1.0 / 10) - (11 * l*l) / 210,0,0,0,-(6.0 / (5 * l)) + (9 * l) / 70,0,-1 + (3 * l) / 5 + (13 * l*l) / 420, },
					{ 0,0,0,(l / 3),0,0,0,0,0,(l / 6),0,0 },
					{ 0,0,-(1.0 / 10) - (11 * l*l) / 210,0,(2 * l) / 15 + l*l*l / 105,0,0,0,1.0 / 10 - (13 * l*l) / 420,0,l / 6 - (3 * l*l) / 20 - l*l*l / 140,0 },
					{ 0,1.0 / 10 + l / 30 + l*l / 30,0,0,0,(2 * l) / 15 + l*l / 105,0,-(1.0 / 10) + (13 * l*l) / 420,0,0,0,l / 6 - (3 * l*l) / 20 - l*l*l / 140 },
					{ (l / 6),0,0,0,0,0,(l / 3),0,0,0,0,0 },
					{ 0,8.0 / 15 - 6.0 / (5 * l) - (3 * l) / 10,0,0,0,-(1.0 / 10) + (13 * l*l) / 420,0,6.0 / (5 * l) + (13 * l) / 35,0,0,0,-1 + (3 * l) / 5 - (11 * l*l) / 210 },
					{ 0,0,-(6.0 / (5 * l)) + (9 * l) / 70,0,1. / 10 - (13 * l*l) / 420,0,0,0,6.0 / (5 * l) + (13 * l) / 35,0,1 - (3 * l) / 5 + (11 * l*l) / 210,0 },
					{ 0,0,0,(l / 6),0,0,0,0,0,(l / 3),0,0 },
					{ 0,0,-1 + (3 * l) / 5 + (13 * l*l) / 420,0,l / 6 - (3 * l*l) / 20 - l*l*l / 140,0,0,0,1 - (3 * l) / 5 + (11 * l*l) / 210,0,(4 * l) / 3 - (12 * l*l) / 5 + (136 * l*l*l) / 105,0 },
					{ 0,1 - (2 * l) / 3 + l*l / 60,0,0,0,l / 6 - (3 * l*l) / 20 - l*l*l / 140,0,-1 + (3 * l) / 5 - (11 * l*l) / 210,0,0,0,(4 * l) / 3 - (12 * l*l) / 5 + (136 * l*l*l) / 105 } };
					//this->damp=Matrix<T>(12,12);
					this->damp = 766 * this->damp;
				}
				void FrictionCoefficientDamp()
				{

				}
				void RayleighDamp()
				{
					//don't need element damp
				}

				void CalculateTransformation()
				{
					T ca = std::cos(this->wellAngle);
					T sa = std::sin(this->wellAngle);
					T cb = std::cos(this->azimuthAngle);
					T sb = std::sin(this->azimuthAngle);
					this->transformation = { {ca,sa*cb,sa*sb,0,0,0,0,0,0,0,0,0},
					{-sa,ca*cb,ca*sb,0,0,0,0,0,0,0,0,0 },
					{0,-sb,cb,0,0,0,0,0,0,0,0,0 },
					{0,0,0,ca,sa*cb,sa*sb ,0,0,0,0,0,0},
					{0,0,0,-sa,ca*cb,ca*sb ,0,0,0,0,0,0},
					{0,0,0,0,-sb,cb ,0,0,0,0,0,0},
					{ 0,0,0,0,0,0,ca,sa*cb,sa*sb ,0,0,0 },
					{ 0,0,0,0,0,0,-sa,ca*cb,ca*sb ,0,0,0 },
					{ 0,0,0,0,0,0,0,-sb,cb ,0,0,0 },
					{ 0,0,0,0,0,0,0,0,0 ,ca,sa*cb,sa*sb },
					{ 0,0,0,0,0,0,0,0,0 ,-sa,ca*cb,ca*sb },
					{ 0,0,0,0,0,0,0,0,0 ,0,-sb,cb } };

				}
			protected:

			};

			template<typename T>
			class SpatialBeamSolver //:public IDynamicSolver<T>
			{
			public:
				SpatialBeamSolver()
				{

				}
				void Solve(std::vector<SpatialBeamElement<T>*>& _element, T _timeInterval, DampType _type = DampType::Rayleigh,
					T _omega1 = 0.2507,T _omega2 = 0.4413,T _zeta1 = 0.2,T _zeta2 = 0.2)
				{
					this->NodeIterationBottomFunctions;

					this->omega1 = _omega1;
					this->omega2 = _omega2;
					this->zeta1 = _zeta1;
					this->zeta2 = _zeta2;

					this->type = _type;
					this->timeInterval = _timeInterval;
					//assemb integral matrix
					this->element = _element;
					int length = this->Length();
					std::vector<int> unknow;
					for (int i = 0; i<length; i++)
					{
						unknow.push_back(0);
					}

					this->mass = Matrix<T>(length, length);
					this->damp = Matrix<T>(length, length);
					this->stiffness = Matrix<T>(length, length);
					for (int i = 0; i < this->element[0]->Force().size(); i++)
					{
						this->force.push_back(Matrix<T>(length, 1));
						this->displacement.push_back(Matrix<T>(length, 1));
						this->velocity.push_back(Matrix<T>(length, 1));
						this->acceleration.push_back(Matrix<T>(length, 1));

						this->unKnowForce.push_back(unknow);
						this->unKnowDisplacement.push_back(unknow);
						this->unKnowVelocity.push_back(unknow);
						this->unKnowAcceleration.push_back(unknow);
					}

					for (IDynamicElement<T>* e : this->element)
					{
						int ei = e->Node()[0];
						int ej = e->Node()[1];

						int ni = ei, nj = ei;
						for (int i = 0; i < e->Stiffness().Row(); i++)
						{
							if (i >= 6)
							{
								ni = ej;
							}
							else
							{
								ni = ei;
							}

							for (int k = 0; k < this->force.size(); k++)
							{
								this->force[k][(ni - 1) * 6 + i % 6][0] += e->Force()[k][i][0];

								this->displacement[k][(ni - 1) * 6 + i % 6][0] += e->Displacement()[k][i][0];
								this->velocity[k][(ni - 1) * 6 + i % 6][0] += e->Velocity()[k][i][0];
								this->acceleration[k][(ni - 1) * 6 + i % 6][0] += e->Acceleration()[k][i][0];

								this->unKnowForce[k][(ni - 1) * 6 + i % 6] += e->UnKnowForce()[k][i];

								this->unKnowDisplacement[k][(ni - 1) * 6 + i % 6] += e->UnKnowDisplacement()[k][i];
								this->unKnowVelocity[k][(ni - 1) * 6 + i % 6] += e->UnKnowVelocity()[k][i];
								this->unKnowAcceleration[k][(ni - 1) * 6 + i % 6] += e->UnKnowAcceleration()[k][i];

							}

							for (int j = 0; j < e->Stiffness().Column(); j++)
							{
								if (j >= 6)
								{
									nj = ej;
								}
								else
								{
									nj = ei;
								}

								this->mass[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Mass()[i][j];
								this->stiffness[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Stiffness()[i][j];
								if (this->type == DampType::Liner)
								{
									this->damp[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Damp()[i][j];
								}
							}
						}
					}
					if (this->type == DampType::Rayleigh)
					{
						this->CalculateDamp(this->type);
					}

					//time--->>>
					this->timeNumber = 0;//bar laji ...
					for (int t = 0; t < this->force.size() - 1; t++)
					{
						if (this->stopSolve == true)
						{
							return;
						}
						this->timeNumber++;//bar laji ...
						this->Newmark(this->stiffness, this->mass, this->damp,
							this->displacement[t], this->velocity[t], this->acceleration[t], this->force[t], this->displacement[t + 1], this->velocity[t + 1], this->acceleration[t + 1], this->force[t + 1],
							this->unKnowDisplacement[t], this->unKnowVelocity[t], this->unKnowAcceleration[t], this->unKnowDisplacement[t + 1], this->unKnowVelocity[t + 1], this->unKnowAcceleration[t + 1], this->unKnowForce[t + 1]);
						
						WriteMatrixToCSV("afterNewmarkV.csv",this->velocity[t+ 1]);
						this->stiffness = Matrix<T>(length, length);
						for (IDynamicElement<T>* e : this->element)
						{
							int ei = e->Node()[0];
							int ej = e->Node()[1];

							int ni = ei, nj = ei;
							for (int i = 0; i < e->Stiffness().Row(); i++)
							{
								if (i >= 6)
								{
									ni = ej;
								}
								else
								{
									ni = ei;
								}

								e->Force()[t][i][0] = this->force[t][(ni - 1) * 6 + i % 6][0];
								e->Displacement()[t][i][0] = this->displacement[t][(ni - 1) * 6 + i % 6][0];
								e->Velocity()[t][i][0] = this->velocity[t][(ni - 1) * 6 + i % 6][0];
								e->Acceleration()[t][i][0] = this->acceleration[t][(ni - 1) * 6 + i % 6][0];

								e->UnKnowForce()[t][i] = this->unKnowForce[t][i];
								e->UnKnowDisplacement()[t][i] = this->unKnowDisplacement[t][i];
								e->UnKnowVelocity()[t][i] = this->unKnowVelocity[t][i];
								e->UnKnowDisplacement()[t][i] = this->unKnowDisplacement[t][i];
							}
							//ReCalculate Stiffness nonliner
							e->DynamicReset(this->type,t);
							WriteMatrixToCSV("eKtemp.csv", e->Stiffness());
							//full Stiffness
							ni = ei, nj = ei;
							for (int i = 0; i < e->Stiffness().Row(); i++)
							{
								if (i >= 6)
								{
									ni = ej;
								}
								else
								{
									ni = ei;
								}

								/*
								for (int k = 0; k < this->force.size(); k++)
								{
									this->force[k][(ni - 1) * 6 + i % 6][0] += e->Force()[k][i][0];

									this->displacement[k][(ni - 1) * 6 + i % 6][0] += e->Displacement()[k][i][0];
									this->velocity[k][(ni - 1) * 6 + i % 6][0] += e->Velocity()[k][i][0];
									this->acceleration[k][(ni - 1) * 6 + i % 6][0] += e->Acceleration()[k][i][0];

									this->unKnowForce[k][(ni - 1) * 6 + i % 6] += e->UnKnowForce()[k][i];

									this->unKnowDisplacement[k][(ni - 1) * 6 + i % 6] += e->UnKnowDisplacement()[k][i];
									this->unKnowVelocity[k][(ni - 1) * 6 + i % 6] += e->UnKnowVelocity()[k][i];
									this->unKnowAcceleration[k][(ni - 1) * 6 + i % 6] += e->UnKnowAcceleration()[k][i];

								}
								*/
								for (int j = 0; j < e->Stiffness().Column(); j++)
								{
									if (j >= 6)
									{
										nj = ej;
									}
									else
									{
										nj = ei;
									}

									this->stiffness[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Stiffness()[i][j];
									if (this->type == DampType::Liner)
									{
										this->damp[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Damp()[i][j];
									}
								}
							}
						}
						if (this->type == DampType::Rayleigh)
						{
							//this->CalculateDamp(this->type);
						}
					}
				}

				std::vector<std::function<void(Matrix<T> &K,Matrix<T> &F,Matrix<T> &d)>> NodeIterationBottomFunctions;


				std::vector<SpatialBeamElement<T>*>& Element()
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

				std::vector<std::vector<int>>& UnKnowForce()
				{
					return this->unKnowForce;
				}
				std::vector<std::vector<int>>& UnKnowDisplacement()
				{
					return this->unKnowDisplacement;
				}
				std::vector<std::vector<int>>& UnKnowVelocity()
				{
					return this->unKnowVelocity;
				}
				std::vector<std::vector<int>>& UnKnowAcceleration()
				{
					return this->unKnowAcceleration;
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

				T TimeInterval()
				{
					return this->timeInterval;
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

				int timeNumber = 0;
				bool stopSolve = false;
			private:
				std::vector<SpatialBeamElement<T>*> element;
				std::vector<Matrix<T>> force;

				std::vector<Matrix<T>> displacement;
				std::vector<Matrix<T>> velocity;
				std::vector<Matrix<T>> acceleration;

				std::vector<std::vector<int>> unKnowForce;
				std::vector<std::vector<int>> unKnowDisplacement;
				std::vector<std::vector<int>> unKnowVelocity;
				std::vector<std::vector<int>> unKnowAcceleration;

				Matrix<T> mass;
				Matrix<T> damp;
				Matrix<T> stiffness;

				T timeInterval;

				DampType type;

				T omega1 = 0.2507;
				T omega2 = 0.4413;
				T zeta1 = 0.2;
				T zeta2 = 0.2;

				/// <summary>
				/// Nodes the iteration.
				/// kd=f
				/// </summary>
				/// <param name="K">The k.</param>
				/// <param name="d">The d.</param>
				/// <param name="F">The f.</param>
				void NodeIteration(Matrix<T>& K, Matrix<T>& d, Matrix<T>& F)
				{
					//3*6=18宽？
					//Three diagonal matrix??
					Matrix<T> d0 = { { d[0][0] },{ d[1][0] },{ d[2][0] },{ d[3][0] },{ d[4][0] },{ d[5][0] } };
					Matrix<T> d1 = { { d[6][0] },{ d[7][0] },{ d[8][0] },{ d[9][0] },{ d[10][0] },{ d[11][0] } };
					Matrix<T> d2 = { { d[12][0] },{ d[13][0] },{ d[14][0] },{ d[15][0] },{ d[16][0] },{ d[17][0] } };

					Matrix<T> F0 = { { F[0][0] },{ F[1][0] },{ F[2][0] },{ F[3][0] },{ F[4][0] },{ F[5][0] } };
					Matrix<T> F1 = { { F[6][0] },{ F[7][0] },{ F[8][0] },{ F[9][0] },{ F[10][0] },{ F[11][0] } };
					Matrix<T> F2 = { { F[12][0] },{ F[13][0] },{ F[14][0] },{ F[15][0] },{ F[16][0] },{ F[17][0] } };

					Matrix<double> K00(6, 6);
					Matrix<double> K01(6, 6);
					Matrix<double> K10(6, 6);
					Matrix<double> K11(6, 6);
					Matrix<double> K12(6, 6);
					Matrix<double> K21(6, 6);
					Matrix<double> K22(6, 6);
					for (int i = 0; i < 6; i++)
					{
						for (int j = 0; j < 6; j++)
						{
							K00[i][j] = K[i][j];
							K01[i][j] = K[i][j + 6];
							K10[i][j] = K[i + 6][j];
							K11[i][j] = K[i + 6][j + 6];
							K12[i][j] = K[i + 6][j + 12];
							K21[i][j] = K[i + 12][j + 6];
							K22[i][j] = K[i + 12][j + 12];
						}
					}

					d1 = (K01^-1)*(F0 - K00*d0);
					for (int i = 0; i < 6; i++)
					{
						d[i + 6][0] = d1[i][0];
					}

					for (int t = 2; t < d.Row()/6; t++)
					{
						for (int j = 0; j < 6; j++)
						{
							d0[j][0] = d[(t - 2) * 6 + j][0];
							d1[j][0] = d[(t - 1) * 6 + j][0];
							d2[j][0] = d[t * 6 + j][0];
							F0[j][0] = F[(t - 2) * 6 + j][0];
							F1[j][0] = F[(t - 1) * 6 + j][0];
							F2[j][0] = F[t * 6 + j][0];
						}

						for (int i = 0; i < 6; i++)
						{
							for (int j = 0; j < 6; j++)
							{
								K00[i][j] = K[(t - 2) * 6 + i][(t - 2) * 6 + j];
								K01[i][j] = K[(t - 2) * 6 + i][(t - 1) * 6 + j];
								K10[i][j] = K[(t - 1) * 6 + i][(t - 2) * 6 + j];
								K11[i][j] = K[(t - 1) * 6 + i][(t - 1) * 6 + j];
								K12[i][j] = K[(t - 1) * 6 + i][t * 6 + j];
								K21[i][j] = K[t * 6 + i][(t - 1) * 6 + j];
								K22[i][j] = K[t * 6 + i][t * 6 + j];
							}
						}

						d2 = (K12^-1)*(F1 - K10*d0 - K11*d1);
						for (int j = 0; j < 6; j++)
						{
							d[t * 6 + j][0] = d2[j][0];
						}

					}
				}
				void NodeIterationBottom(Matrix<T>& K, Matrix<T>& d, Matrix<T>& F)
				{

					int r = d.Row();
					//3*6=18宽？
					//Three diagonal matrix??
					Matrix<T> d0 = { { d[r - 18][0] },{ d[r - 17][0] },{ d[r - 16][0] },{ d[r - 15][0] },{ d[r - 14][0] },{ d[r - 13][0] } };
					Matrix<T> d1 = { { d[r - 12][0] },{ d[r - 11][0] },{ d[r - 10][0] },{ d[r - 9][0] },{ d[r - 8][0] },{ d[r - 7][0] } };
					Matrix<T> d2 = { { d[r - 6][0] },{ d[r - 5][0] },{ d[r - 4][0] },{ d[r - 3][0] },{ d[r - 2][0] },{ d[r - 1][0] } };

					Matrix<T> F0 = { { F[r - 18][0] },{ F[r - 17][0] },{ F[r - 16][0] },{ F[r - 15][0] },{ F[r - 14][0] },{ F[r - 13][0] } };
					Matrix<T> F1 = { { F[r - 12][0] },{ F[r - 11][0] },{ F[r - 10][0] },{ F[r - 9][0] },{ F[r - 8][0] },{ F[r - 7][0] } };
					Matrix<T> F2 = { { F[r - 6][0] },{ F[r - 5][0] },{ F[r - 4][0] },{ F[r - 3][0] },{ F[r - 2][0] },{ F[r - 1][0] } };

					Matrix<double> K00(6, 6);
					Matrix<double> K01(6, 6);
					Matrix<double> K10(6, 6);
					Matrix<double> K11(6, 6);
					Matrix<double> K12(6, 6);
					Matrix<double> K21(6, 6);
					Matrix<double> K22(6, 6);
					for (int i = 0; i < 6; i++)
					{
						for (int j = 0; j < 6; j++)
						{
							K00[i][j] = K[i + r - 18][j + r - 18];
							K01[i][j] = K[i + r - 18][j + r - 12];
							K10[i][j] = K[i + r - 12][j + r - 18];
							K11[i][j] = K[i + r - 12][j + r - 12];
							K12[i][j] = K[i + r - 12][j + r - 6];
							K21[i][j] = K[i + r - 6][j + r - 12];
							K22[i][j] = K[i + r - 6][j + r - 6];
						}
					}

					d1 = (K21^-1)*(F2 - K22*d2);
					for (int i = 0; i < 6; i++)
					{
						d[i + r - 12][0] = d1[i][0];
					}
					//Top force
					for (auto df : this->NodeIterationBottomFunctions)
					{
						df(K, F, d);
					}
					//Top force
					/*Matrix<T> Ftop = K21*d1 + K22*d2;
					for (int i = 0; i < Ftop.Row(); i++)
					{
						F[i][0] = Ftop[i][0];
					}*/

					for (int t = r / 6 - 1; t >= 2; t--)
					{
						for (int j = 0; j < 6; j++)
						{
							d0[j][0] = d[(t - 2) * 6 + j][0];
							d1[j][0] = d[(t - 1) * 6 + j][0];
							d2[j][0] = d[t * 6 + j][0];
							F0[j][0] = F[(t - 2) * 6 + j][0];
							F1[j][0] = F[(t - 1) * 6 + j][0];
							F2[j][0] = F[t * 6 + j][0];
						}

						for (int i = 0; i < 6; i++)
						{
							for (int j = 0; j < 6; j++)
							{
								K00[i][j] = K[(t - 2) * 6 + i][(t - 2) * 6 + j];
								K01[i][j] = K[(t - 2) * 6 + i][(t - 1) * 6 + j];
								K10[i][j] = K[(t - 1) * 6 + i][(t - 2) * 6 + j];
								K11[i][j] = K[(t - 1) * 6 + i][(t - 1) * 6 + j];
								K12[i][j] = K[(t - 1) * 6 + i][t * 6 + j];
								K21[i][j] = K[t * 6 + i][(t - 1) * 6 + j];
								K22[i][j] = K[t * 6 + i][t * 6 + j];
							}
						}

						d0 = (K10^-1)*(F1 - K11*d1 - K12*d2);
						for (int j = 0; j < 6; j++)
						{
							d[(t - 2) * 6 + j][0] = d0[j][0];
						}

					}
				}
				void Newmark(Matrix<T>& K, Matrix<T>& M, Matrix<T>& C,
					Matrix<T>& d, Matrix<T>& v, Matrix<T>& a,Matrix<T>& F, Matrix<T>& dnext, Matrix<T>& vnext, Matrix<T>& anext, Matrix<T>& Fnext,
					std::vector<int>& unKnowD, std::vector<int>& unKnowV, std::vector<int>& unKnowA, std::vector<int>& unKnowDnext, std::vector<int>& unKnowVnext, std::vector<int>& unKnowAnext, std::vector<int>& unKnowFnext)
				{
					T DELTAt = this->timeInterval;
					T eta = 0.25;
					T delta = 0.5;
					T c[8] = { 1.0 / eta / DELTAt / DELTAt ,delta / eta / DELTAt ,1.0 / eta / DELTAt ,1.0 / 2 / eta - 1,delta / eta - 1, DELTAt*(delta / 2 / eta - 1), DELTAt*(1 - delta),delta*DELTAt };
					Matrix<T> effectiveK = K + c[0] * M + c[1] * C;
					Matrix<T> effectiveFnext = Fnext + M*(c[0] * d + c[2] * v + c[4] * a) + C*(c[1] * d + c[3] * v + c[5] * a);

					//d(t+DELTAt),v(t+DELTAt),a(t+DELTAt)
					//this->NodeIteration(effectiveK, dnext, effectiveFnext);

					//CutZeroAddOne(effectiveK,dnext,effectiveFnext,unKnowDnext,unKnowFnext,this->Length());
					//WriteMatrixToCSV("M.csv", M);
					WriteMatrixToCSV("eK.csv", effectiveK);
					WriteMatrixToCSV("eFn.csv", effectiveFnext);
					//WriteMatrixToCSV("eKi.csv", effectiveK.Inverse());
					//dnext = effectiveK.Inverse()*effectiveFnext;
					auto nCut = [&]() {
						int r = dnext.Row();
						//3*6=18宽？
						//Three diagonal matrix??
						Matrix<T> d0 = { { dnext[r - 18][0] },{ dnext[r - 17][0] },{ dnext[r - 16][0] },{ dnext[r - 15][0] },{ dnext[r - 14][0] },{ dnext[r - 13][0] } };
						Matrix<T> d1 = { { dnext[r - 12][0] },{ dnext[r - 11][0] },{ dnext[r - 10][0] },{ dnext[r - 9][0] },{ dnext[r - 8][0] },{ dnext[r - 7][0] } };
						Matrix<T> d2 = { { dnext[r - 6][0] },{ dnext[r - 5][0] },{ dnext[r - 4][0] },{ dnext[r - 3][0] },{ dnext[r - 2][0] },{ dnext[r - 1][0] } };

						Matrix<T> F0 = { { effectiveFnext[r - 18][0] },{ effectiveFnext[r - 17][0] },{ effectiveFnext[r - 16][0] },{ effectiveFnext[r - 15][0] },{ effectiveFnext[r - 14][0] },{ effectiveFnext[r - 13][0] } };
						Matrix<T> F1 = { { effectiveFnext[r - 12][0] },{ effectiveFnext[r - 11][0] },{ effectiveFnext[r - 10][0] },{ effectiveFnext[r - 9][0] },{ effectiveFnext[r - 8][0] },{ effectiveFnext[r - 7][0] } };
						Matrix<T> F2 = { { effectiveFnext[r - 6][0] },{ effectiveFnext[r - 5][0] },{ effectiveFnext[r - 4][0] },{ effectiveFnext[r - 3][0] },{ effectiveFnext[r - 2][0] },{ effectiveFnext[r - 1][0] } };

						Matrix<double> K00(6, 6);
						Matrix<double> K01(6, 6);
						Matrix<double> K10(6, 6);
						Matrix<double> K11(6, 6);
						Matrix<double> K12(6, 6);
						Matrix<double> K21(6, 6);
						Matrix<double> K22(6, 6);
						for (int i = 0; i < 6; i++)
						{
							for (int j = 0; j < 6; j++)
							{
								K00[i][j] = effectiveK[i + r - 18][j + r - 18];
								K01[i][j] = effectiveK[i + r - 18][j + r - 12];
								K10[i][j] = effectiveK[i + r - 12][j + r - 18];
								K11[i][j] = effectiveK[i + r - 12][j + r - 12];
								K12[i][j] = effectiveK[i + r - 12][j + r - 6];
								K21[i][j] = effectiveK[i + r - 6][j + r - 12];
								K22[i][j] = effectiveK[i + r - 6][j + r - 6];
							}
						}

						d1 = (K21^-1)*(F2 - K22*d2);
						for (int i = 0; i < 6; i++)
						{
							dnext[i + r - 12][0] = d1[i][0];
						}
						for (auto df : this->NodeIterationBottomFunctions)
						{
							df(K, Fnext, dnext);
						}
						effectiveFnext = Fnext + M*(c[0] * d + c[2] * v + c[4] * a) + C*(c[1] * d + c[3] * v + c[5] * a);

						Matrix<T> eKz(effectiveK.Row() - 6, effectiveK.Column() - 6);
						Matrix<T> eFnz(effectiveFnext.Row() - 6, effectiveFnext.Column());

						for (int i = 0; i < eKz.Row(); i++)
						{
							eFnz[i][0] = effectiveFnext[i][0];
							for (int j = 0; j < eKz.Column(); j++)
							{
								eKz[i][j] = effectiveK[i][j];
							}
						}

						Matrix<T> dnz = (eKz^-1)*eFnz;

						for (int i = 0; i < dnz.Row(); i++)
						{
							dnext[i][0] = dnz[i][0];
						}
					};
					auto nBottom = [&]()->void {
						int r = dnext.Row();
						//3*6=18宽？
						//Three diagonal matrix??
						Matrix<T> d0 = { { dnext[r - 18][0] },{ dnext[r - 17][0] },{ dnext[r - 16][0] },{ dnext[r - 15][0] },{ dnext[r - 14][0] },{ dnext[r - 13][0] } };
						Matrix<T> d1 = { { dnext[r - 12][0] },{ dnext[r - 11][0] },{ dnext[r - 10][0] },{ dnext[r - 9][0] },{ dnext[r - 8][0] },{ dnext[r - 7][0] } };
						Matrix<T> d2 = { { dnext[r - 6][0] },{ dnext[r - 5][0] },{ dnext[r - 4][0] },{ dnext[r - 3][0] },{ dnext[r - 2][0] },{ dnext[r - 1][0] } };

						Matrix<T> F0 = { { effectiveFnext[r - 18][0] },{ effectiveFnext[r - 17][0] },{ effectiveFnext[r - 16][0] },{ effectiveFnext[r - 15][0] },{ effectiveFnext[r - 14][0] },{ effectiveFnext[r - 13][0] } };
						Matrix<T> F1 = { { effectiveFnext[r - 12][0] },{ effectiveFnext[r - 11][0] },{ effectiveFnext[r - 10][0] },{ effectiveFnext[r - 9][0] },{ effectiveFnext[r - 8][0] },{ effectiveFnext[r - 7][0] } };
						Matrix<T> F2 = { { effectiveFnext[r - 6][0] },{ effectiveFnext[r - 5][0] },{ effectiveFnext[r - 4][0] },{ effectiveFnext[r - 3][0] },{ effectiveFnext[r - 2][0] },{ effectiveFnext[r - 1][0] } };

						Matrix<double> K00(6, 6);
						Matrix<double> K01(6, 6);
						Matrix<double> K10(6, 6);
						Matrix<double> K11(6, 6);
						Matrix<double> K12(6, 6);
						Matrix<double> K21(6, 6);
						Matrix<double> K22(6, 6);
						for (int i = 0; i < 6; i++)
						{
							for (int j = 0; j < 6; j++)
							{
								K00[i][j] = effectiveK[i + r - 18][j + r - 18];
								K01[i][j] = effectiveK[i + r - 18][j + r - 12];
								K10[i][j] = effectiveK[i + r - 12][j + r - 18];
								K11[i][j] = effectiveK[i + r - 12][j + r - 12];
								K12[i][j] = effectiveK[i + r - 12][j + r - 6];
								K21[i][j] = effectiveK[i + r - 6][j + r - 12];
								K22[i][j] = effectiveK[i + r - 6][j + r - 6];
							}
						}

						d1 = (K21^-1)*(F2 - K22*d2);
						for (int i = 0; i < 6; i++)
						{
							dnext[i + r - 12][0] = d1[i][0];
						}
						//Top force
						for (auto df : this->NodeIterationBottomFunctions)
						{
							df(K, Fnext, dnext);
						}
						effectiveFnext = Fnext + M*(c[0] * d + c[2] * v + c[4] * a) + C*(c[1] * d + c[3] * v + c[5] * a);
						F0 = { { effectiveFnext[r - 18][0] },{ effectiveFnext[r - 17][0] },{ effectiveFnext[r - 16][0] },{ effectiveFnext[r - 15][0] },{ effectiveFnext[r - 14][0] },{ effectiveFnext[r - 13][0] } };
						F1 = { { effectiveFnext[r - 12][0] },{ effectiveFnext[r - 11][0] },{ effectiveFnext[r - 10][0] },{ effectiveFnext[r - 9][0] },{ effectiveFnext[r - 8][0] },{ effectiveFnext[r - 7][0] } };
						F2 = { { effectiveFnext[r - 6][0] },{ effectiveFnext[r - 5][0] },{ effectiveFnext[r - 4][0] },{ effectiveFnext[r - 3][0] },{ effectiveFnext[r - 2][0] },{ effectiveFnext[r - 1][0] } };
						//Top force
						/*Matrix<T> Ftop = K21*d1 + K22*d2;
						for (int i = 0; i < Ftop.Row(); i++)
						{
						effectiveFnext[i][0] = Ftop[i][0];
						}*/

						Matrix<T> cuteK(r - 2 * 6, r - 2 * 6);
						Matrix<T> cuteFnext(r - 2 * 6, 1);
						for (int i = 0; i < cuteK.Row(); i++)
						{
							cuteFnext[i][0] = effectiveFnext[i][0];
							for (int j = 0; j < cuteK.Column(); j++)
							{
								cuteK[i][j] = effectiveK[i][j];
							}
						}

						Matrix<T> cutednext = (cuteK^-1)*cuteFnext;
						for (int i = 0; i < cutednext.Row(); i++)
						{
							dnext[i][0] = cutednext[i][0];
						}

						/*

						for (int t = r / 6 - 1; t >= 2; t--)
						{
							for (int j = 0; j < 6; j++)
							{
								d0[j][0] = dnext[(t - 2) * 6 + j][0];
								d1[j][0] = dnext[(t - 1) * 6 + j][0];
								d2[j][0] = dnext[t * 6 + j][0];
								F0[j][0] = effectiveFnext[(t - 2) * 6 + j][0];
								F1[j][0] = effectiveFnext[(t - 1) * 6 + j][0];
								F2[j][0] = effectiveFnext[t * 6 + j][0];
							}

							for (int i = 0; i < 6; i++)
							{
								for (int j = 0; j < 6; j++)
								{
									K00[i][j] = effectiveK[(t - 2) * 6 + i][(t - 2) * 6 + j];
									K01[i][j] = effectiveK[(t - 2) * 6 + i][(t - 1) * 6 + j];
									K10[i][j] = effectiveK[(t - 1) * 6 + i][(t - 2) * 6 + j];
									K11[i][j] = effectiveK[(t - 1) * 6 + i][(t - 1) * 6 + j];
									K12[i][j] = effectiveK[(t - 1) * 6 + i][t * 6 + j];
									K21[i][j] = effectiveK[t * 6 + i][(t - 1) * 6 + j];
									K22[i][j] = effectiveK[t * 6 + i][t * 6 + j];
								}
							}

							d0 = (K10^-1)*(F1 - K11*d1 - K12*d2);
							for (int j = 0; j < 6; j++)
							{
								dnext[(t - 2) * 6 + j][0] = d0[j][0];
							}
							
							if (t == 2)
							{
								d0 = (K00^-1)*(F0 - K01*d1);
								for (int j = 0; j < 6; j++)
								{
									dnext[(t - 2) * 6 + j][0] = d0[j][0];
								}
							}
						}
						*/					
					};
					//nBottom();
					//nCut();
					//this->NodeIterationBottom(effectiveK, dnext, effectiveFnext);
					this->NodeIteration(effectiveK, dnext, effectiveFnext);

					auto nFast = [&]() {
						anext = c[0] * (dnext - d) - c[2] * v - c[3] * a;
						vnext = v + c[6] * a + c[7] * anext;
					};

					auto vFast = [&]()->void {
						Matrix<T> effectiveKv = 1.0 / delta / DELTAt*M + C;
						Matrix<T> effectiveFnextv = Fnext + M*(1.0 / delta / DELTAt*v + (1.0 / delta + 1)*a) - K*dnext;

						WriteMatrixToCSV("effKv.csv", effectiveKv);
						WriteMatrixToCSV("effFv.csv", effectiveFnextv);

						this->NodeIteration(effectiveKv, vnext, effectiveFnextv);

						//WriteMatrixToCSV("effVnext.csv", vnext);

						Matrix<T> effectiveKa = M;
						Matrix<T> effectiveFnexta = Fnext - C*vnext - K*dnext;

						WriteMatrixToCSV("effKa.csv", effectiveKa);
						WriteMatrixToCSV("effFa.csv", effectiveFnexta);

						this->NodeIteration(effectiveKa, anext, effectiveFnexta);
					};
					auto aFast = [&]()->void {
						Matrix<T> effectiveKa = M + DELTAt*delta*C;
						Matrix<T> effectiveFnexta = Fnext - C*v - DELTAt*(1 + delta)*C*a - K*dnext;

						WriteMatrixToCSV("effKa.csv", effectiveKa);
						WriteMatrixToCSV("effFa.csv", effectiveFnexta);

						this->NodeIteration(effectiveKa, anext, effectiveFnexta);

						//WriteMatrixToCSV("effVnext.csv", vnext);

						Matrix<T> effectiveKv = C;
						Matrix<T> effectiveFnextv = Fnext - M*anext - K*dnext;

						WriteMatrixToCSV("effKv.csv", effectiveKv);
						WriteMatrixToCSV("effFv.csv", effectiveFnextv);

						this->NodeIteration(effectiveKv, vnext, effectiveFnextv);
					};
					
					nFast();
					/*
					if (this->ifFirstNewMark == true)
					{
						v = 1 / DELTAt*(dnext - d);
						a = (M^-1) * (F - C*v - K*d);
					}
					*/

					//WriteMatrixToCSV("anNewmarkDnext.csv", dnext);
					//temp=c[0] * (dnext - d) - c[2] * v - c[3] * a;
					//CutZeroAddOne(EM,anext,temp,unKnowAnext,std::vector<int>(this->Length(),1),this->Length());
					//
					//temp=v + ((1 - delta)*a + delta*anext)*DELTAt;
					//CutZeroAddOne(EM,vnext,temp,unKnowVnext,std::vector<int>(this->Length(),1),this->Length());
					//
					//temp=M*anext+C*vnext+K*dnext;
					//CutZeroAddOne(EM,Fnext,temp,unKnowFnext,std::vector<int>(this->Length(),1),this->Length());
					//Fnext=M*anext+C*vnext+K*dnext;

					this->ifFirstNewMark = false;
				}

				void CalculateDamp(DampType type)
				{
					switch (type) {
					case DampType::Liner:
						this->LinerDamp();
						break;
					case DampType::FrictionCoefficient:
						this->FrictionCoefficientDamp();
						break;
					case DampType::Rayleigh:
						this->RayleighDamp();
						break;
					default:
						this->RayleighDamp();
						break;
					}
				}
				void LinerDamp()
				{

				}
				void FrictionCoefficientDamp()
				{
					T omega1 = 1, omega2 = 2;
					T zeta1 = 0.01, zeta2 = 0.02;

					T phi = 2 * (zeta1*omega2 - zeta2*omega1) / (omega2*omega2 - omega1*omega1);
					T xi = omega1*omega2*phi;

					this->damp = xi*this->mass + phi*this->stiffness;
				}
				void RayleighDamp()
				{
					//WriteMatrixToCSV("K.csv",this->stiffness);
					//WriteMatrixToCSV("M.csv",this->mass);
					//Matrix<T> p=JacobiNaturalFrequency(this->stiffness,this->mass);
					//WriteMatrixToCSV("p.csv",p);
					//T omega1=p[0][0], omega2=p[1][0];

					T phi = 2 * (zeta1*omega2 - zeta2*omega1) / (omega2*omega2 - omega1*omega1);
					T xi = omega1*omega2*phi;

					this->damp = xi*this->mass + phi*this->stiffness;
				}

				bool ifFirstNewMark = true;
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

			template<typename T>
			class SpatialBeamPreprocessor
			{
			public:
				SpatialBeamElement<T> *element;

				SpatialBeamSolver<T> *solver;

				T g = 9.8;
				T densitySuckerRod;
				T densityCentralizer;
				
				T sectionalAreaSuckerRod;
				T sectionalAreaCentralizer;
				T lengthSuckerRod;
				T lengthCentralizer;

				T dynamicViscosity;
				T densityWellFluid;
				T wellDiameter;
				T timeInterval;

				Matrix<T> gravitySuckerRod;

				SpatialBeamPreprocessor(T _wellDiameter,T _densityWellFluid,T _dynamicViscosity,T _timeInterval)
				{
					this->wellDiameter = _wellDiameter;
					this->densityWellFluid = _densityWellFluid;
					this->dynamicViscosity = _dynamicViscosity;
					this->timeInterval = _timeInterval;
				}

				//浮重：Gravity&flo
				void Gravity()
				{
					T q1 = (this->element->Density() - this->densityWellFluid)*this->g*this->element->SectionArea()*std::cos(this->element->WellAngle());
					T q2 = (this->element->Density() - this->densityWellFluid)*this->g*this->element->SectionArea()*std::sin(this->element->WellAngle());
					T l = this->element->Length();

					Matrix<T> q = { { -q1*l / 2 },{ -q2*l / 2 },{ 0 },{ 0 },{ 0 },{ -q2*l*l / 12 },{ -q1*l / 2 },{ -q2*l / 2 },{ 0 },{ 0 },{ 0 },{ q2*l*l / 12 } };
					for (Matrix<T> &F : this->element->Force())
					{
						F = F + q;
					}
				}

				//should use it after Gravity.
				void GravityStatic(std::vector<SpatialBeamElement<T>*>& _element)
				{
					//consecutive node number : 1,2,3,4,...
					int all = 0;
					for (IDynamicElement<T>* e : _element)
					{
						for (int i : e->Node())
						{
							if (i > all)
							{
								all = i;
							}
						}
					}
					int length=all * 6;

					Matrix<T> stiffness(length, length);
					Matrix<T> force(length, 1);
					Matrix<T> displacement(length, 1);

					for (IDynamicElement<T>* e : _element)
					{
						int ei = e->Node()[0];
						int ej = e->Node()[1];

						int ni = ei, nj = ei;
						for (int i = 0; i < e->Stiffness().Row(); i++)
						{
							if (i >= 6)
							{
								ni = ej;
							}
							else
							{
								ni = ei;
							}

							force[(ni - 1) * 6 + i % 6][0] += e->Force()[0][i][0];
							displacement[(ni - 1) * 6 + i % 6][0] += e->Displacement()[0][i][0];

							for (int j = 0; j < e->Stiffness().Column(); j++)
							{
								if (j >= 6)
								{
									nj = ej;
								}
								else
								{
									nj = ei;
								}

								stiffness[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Stiffness()[i][j];
							}
						}
					}

					std::vector<int> unKnowDisplacement;
					for (int i = 0; i < 6; i++)
					{
						unKnowDisplacement.push_back(1);
					}			
					for (int i = 6; i < length; i++)
					{
						unKnowDisplacement.push_back(0);
					}
					std::vector<int> unKnowForce;
					for (int i = 0; i < 6; i++)
					{
						unKnowForce.push_back(0);
					}
					for (int i = 6; i < length; i++)
					{
						unKnowForce.push_back(1);
					}

					//EliminationMethod(stiffness, displacement, force, unKnowDisplacement, unKnowForce, length);
					CutZeroAddOne(stiffness, displacement, force, unKnowDisplacement, unKnowForce, length);

					for (IDynamicElement<T>* e : _element)
					{
						int ei = e->Node()[0];
						int ej = e->Node()[1];

						int ni = ei, nj = ei;
						for (int i = 0; i < e->Stiffness().Row(); i++)
						{
							if (i >= 6)
							{
								ni = ej;
							}
							else
							{
								ni = ei;
							}

							e->Displacement()[0][i][0] = displacement[(ni - 1) * 6 + i % 6][0];
							e->Force()[0][i][0] = force[(ni - 1) * 6 + i % 6][0];
						}
					}

				}

				//液体阻力:Ff
				void LiquidResistance(T At)
				{
					this->element->DynamicResetFunctions.push_back([&](SpatialBeamElement<T> *e,int t) {
						Matrix<T> &v = e->Velocity()[t];
						Matrix<T> &F = e->Force()[t];

						T m = std::sqrt(At / this->element->SectionArea());
						Matrix<T> nv(v.Row(), v.Column());
						for (int i = 0; i<v.Row(); i++)
						{
							nv[i][0] = std::abs(v[i][0]);
						}
						Matrix<T> Ff = -M_PI*this->dynamicViscosity*(m*m - 1) / ((m*m + 1)*std::log(m) / std::log(M_E) - m*m + 1)*(v + nv)*this->element->Length();

						e->Force()[t] = e->Force()[t] + Ff;
					});

				}
				/**
				*示功图加入，按时间排列，两行，第一行力，第二行位移
				*/
				void IndicatorDagram(Matrix<T> &m)
				{
					Matrix<T> Dt(2, m.Column());
					for (int t = 0; t < m.Column(); t++)
					{
						Dt[0][t] = m[1][t];
						Dt[1][t] = t*this->timeInterval;
					}
					Matrix<T> Vt(2, m.Column());
					Vt = InterpolationDerivation(Dt);
					Matrix<T> At(2, m.Column());
					At = InterpolationDerivation(Vt);

					for (int t = 0; t < m.Column(); t++)
					{
						this->element->Force()[t][0][0] = m[0][t];
						this->element->Displacement()[t][0][0] = m[1][t];
						this->element->Velocity()[t][0][0] = Vt[0][t];
						this->element->Acceleration()[t][0][0] = At[0][t];
					}

					//WriteMatrixToCSV("Atop.csv", At);
					
				}
				/**
				*泵功图加入，按时间排列，两行，第一行力，第二行位移
				*/
				void PumpWork(Matrix<T> &m)
				{
					Matrix<T> Dt(2, m.Column());
					for (int t = 0; t < m.Column(); t++)
					{
						Dt[0][t] = m[1][t];
						Dt[1][t] = t*this->timeInterval;
					}
					Matrix<T> Vt(2, m.Column());
					Vt = InterpolationDerivation(Dt);
					Matrix<T> At(2, m.Column());
					At = InterpolationDerivation(Vt);
					for (int t = 0; t < m.Column(); t++)
					{
						//this->element->Force()[t][6][0] = m[0][t];
						this->element->Displacement()[t][6][0] = m[1][t];
						this->element->Velocity()[t][6][0] = Vt[0][t];
						this->element->Acceleration()[t][6][0] = At[0][t];
					}
					int t = 0;
					this->solver->NodeIterationBottomFunctions.push_back([&](Matrix<T> &K, Matrix<T> &F, Matrix<T> &d) {
						int n = F.Row();

						for (int t = 0; t < m.Column(); t++)
						{
							T fn = (this->solver->Stiffness() * this->solver->Displacement()[t])[n - 6][0];

							this->solver->Force()[t][0][0] = m[0][t] - fn;
						}

					});



					/*
					for (int t = 0; t < m.Column(); t++)
					{
						Matrix<T> Fn = this->element->Stiffness()*this->element->Displacement()[t];

						//only index 0?
						top->Force()[t][0][0] = m[0][t] - Fn[6][0];
					}
					*/
					WriteMatrixToCSV("Aend.csv", At);

				}

				void ImpactForce(T kn = 0.008, T ks = 0.0002)
				{
					this->element->DynamicResetFunctions.push_back([&](SpatialBeamElement<T> *e, int t) {
						bool impact = false;
						Matrix<T> k(12, 12);
						k = { { kn,0,0,0,0,0,-kn,0,0,0,0,0 },
						{ 0,ks,0,0,0,0,0,-ks,0,0,0,0 },
						{ 0,0,ks,0,0,0,0,0,-ks,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 },
						{ -kn,0,0,0,0,0,kn,0,0,0,0,0 },
						{ 0,-ks,0,0,0,0,0,ks,0,0,0,0 },
						{ 0,0,-ks,0,0,0,0,0,ks,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 },
						{ 0,0,0,0,0,0,0,0,0,0,0,0 } };
						//直井？
						impact = ((e->Displacement()[t][1][0] * e->Displacement()[t][1][0] + e->Displacement()[t][2][0] * e->Displacement()[t][2][0]) >= (this->wellDiameter*this->wellDiameter / 4)) || ((e->Displacement()[t][7][0] * e->Displacement()[t][7][0] + e->Displacement()[t][8][0] * e->Displacement()[t][8][0]) >= (this->wellDiameter*this->wellDiameter / 4));

						if (impact == true)
						{
							e->Force()[t] = e->Force()[t] + k*(e->Displacement()[t]);
							e->Stiffness() = e->Stiffness() + k;
						}
					});
				}

			};

			template<typename T>
			class SpatialBeamPostprocessor
			{
			public:
				SpatialBeamSolver<T> *solver;

				SpatialBeamPostprocessor(SpatialBeamSolver<T> *_solver)
				{
					this->solver = _solver;
				}

				Matrix<T> PumpWorkBuild()
				{
					WriteMatrixToCSV("Stiffness.csv",this->solver->Stiffness());
					Matrix<T> fs(2, this->solver->Force().size());
					for (int t = 0; t < fs.Column(); t++)
					{
						Matrix<T> Fn = (this->solver->Stiffness())*(this->solver->Displacement()[t]);
						fs[0][t] = this->solver->Force()[t][0][0] + Fn[this->solver->Force()[0].Row() - 6][0];
						fs[1][t] = this->solver->Displacement()[t][this->solver->Force()[0].Row() - 6][0];
					}
					
					return fs;
				}

				Matrix<T> CrankRotatingSpeed(T R1,T R2,T r)
				{
					Matrix<T> rotatingSpeed(2, this->solver->Velocity()->size());

					for (int t = 0; t < this->solver->Velocity()->size(); t++)
					{
						rotatingSpeed[0][t] = t*this->solver->TimeInterval();
						rotatingSpeed[1][t] = 30 * R1*this->solver->Velocity()[t][0][0] / M_PI / R2 / r;
					}
					
					return rotatingSpeed;
				}

				Matrix<T> SuspensionDisplacementToMotorAngle(T r, T R1, T R2, T H, T W, T L, Matrix<T> displacement)
				{
					Matrix<T> angle = displacement;
					for (int i = 0; i < displacement.Column(); i++)
					{
						T d = displacement[1][i];
						T phi0 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L - r)*(L - r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T phi = d / R2 + phi0;
						T H1 = H - R1*std::cos(phi);
						T W1 = W - R1*std::sin(phi);
						T theta1 = std::acos((H1*H1 - L*L + W1*W1 + r*r) / (2 * r*std::sqrt(W1*W1 + H1*H1)));
						T theta2 = std::atan(H1 / W1);

						if(displacement[1][i-1]<d|| i == 0 )
						{
							angle[1][i] = M_PI/2 - theta1 + theta2;
						}
						else
						{
							angle[1][i] = M_PI/2 + theta1 + theta2;
						}

					}

					return angle;
				}

				Matrix<T> Differential(Matrix<T> fx)
				{
					Matrix<T> dfx = fx;
					//back 2 points
					dfx[1][0] = (fx[1][1] - fx[1][0]) / (fx[0][1] - fx[0][0]);
					for (int i = 1; i < fx.Column() - 1; i++)
					{
						//3 points
						dfx[1][i] = (fx[1][i + 1] - fx[1][i - 1]) / (fx[0][i + 1] - fx[0][i - 1]);
					}
					//front 2 points
					dfx[1][fx.Column()-1]= (fx[1][fx.Column() - 1] - fx[1][fx.Column() - 2]) / (fx[0][fx.Column() - 1] - fx[0][fx.Column() - 2]);

					return dfx;
				}

				Multinomial<T> TaylorExpandWalkingBeamAngleToMotorRotationalVelocity(T r, T R1, T R2, T H, T W, T L,bool up=true)
				{
					auto fUp = [&](T phi)->T {
						T phi0 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L - r)*(L - r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T phi1 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L + r)*(L + r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T H1 = H - R1*std::sin(M_PI / 2 - phi);
						T W1 = W - R1*std::cos(M_PI / 2 - phi);
						T theta1 = std::acos((H1*H1 - L*L + W1*W1 + r*r) / (2 * r*std::sqrt(W1*W1 + H1*H1)));
						T theta2 = std::atan(H1 / W1);

						return M_PI / 2 - theta1 + theta2;
					};
					auto fDown = [&](T phi)->T {
						T phi0 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L - r)*(L - r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T phi1 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L + r)*(L + r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T H1 = H + R1*std::sin(phi - M_PI / 2);
						T W1 = W - R1*std::sin(phi - M_PI / 2);
						T theta1 = std::acos((H1*H1 - L*L + W1*W1 + r*r) / (2 * r*std::sqrt(W1*W1 + H1*H1)));
						T theta2 = std::atan(H1 / W1);

						return M_PI / 2 + theta1 + theta2;
					};

					T expandXUp = 1;
					T expandXDown = 1;
					int expandN = 4;
					T deltaX = 0.0001;
					Multinomial<T> result(expandN);
					std::function<T(T, int, std::function<T(T)>, T)> dfn;
					dfn = [&dfn](T x, int n,std::function<T(T)> f,T dx)->T {
						if (n == 0)
						{
							return f(x);
						}
						else
						{
							return (dfn(x + dx, n - 1, f,dx) - dfn(x, n - 1, f,dx)) / dx;
						}
					};

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

					for (int i = 0; i < expandN; i++)
					{
						if (up == true)
						{
							result.coefficients[i] = dfn(expandXUp, i, fUp, deltaX) / nnn(i);
						}
						else
						{
							result.coefficients[i] = dfn(expandXDown, i, fDown, deltaX) / nnn(i);
						}
					}

					return result;
				}
			
				static void MotorRotationalVelocityWithWalkingBeamAngle(T r, T R1, T R2, T H, T W, T L, Matrix<T> phi_t_Up, Matrix<T> phi_t_Down, Matrix<T> &omega_phi_Up, Matrix<T> &omega_phi_Down)
				{
					auto fUp = [&](T phi)->T {
						T phi0 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L - r)*(L - r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T phi1 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L + r)*(L + r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T H1 = H - R1*std::sin(M_PI / 2 - phi);
						T W1 = W - R1*std::cos(M_PI / 2 - phi);
						T theta1 = std::acos((H1*H1 - L*L + W1*W1 + r*r) / (2 * r*std::sqrt(W1*W1 + H1*H1)));
						T theta2 = std::atan(H1 / W1);

						return M_PI / 2 - theta1 + theta2;
					};
					auto fDown = [&](T phi)->T {
						T phi0 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L - r)*(L - r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T phi1 = std::atan(W / H) + std::acos((W*W + H*H + R1*R1 - (L + r)*(L + r)) / (2 * R1*std::sqrt(W*W + H*H)));
						T H1 = H + R1*std::sin(phi - M_PI / 2);
						T W1 = W - R1*std::sin(phi - M_PI / 2);
						T theta1 = std::acos((H1*H1 - L*L + W1*W1 + r*r) / (2 * r*std::sqrt(W1*W1 + H1*H1)));
						T theta2 = std::atan(H1 / W1);

						return M_PI / 2 + theta1 + theta2;
					};


					SpatialBeamSolver<T>* spsolver = nullptr;
					SpatialBeamPostprocessor<T> sppost(spsolver);
					//Multinomial<T> theta_phi_Up = sppost.TaylorExpandWalkingBeamAngleToMotorRotationalVelocity(r, R1, R2, H, W, L);
					//Multinomial<T> theta_phi_Down = sppost.TaylorExpandWalkingBeamAngleToMotorRotationalVelocity(r, R1, R2, H, W, L, false);

					Matrix<T> theta_t_Up = phi_t_Up;
					Matrix<T> theta_t_Down = phi_t_Down;

					for (int i = 0; i < phi_t_Up.Column(); i++)
					{
						theta_t_Up[1][i] = fUp(phi_t_Up[1][i]);
					}
					for (int i = 0; i < phi_t_Down.Column(); i++)
					{
						theta_t_Down[1][i] = fDown(phi_t_Down[1][i]);
					}
					Matrix<T> omega_t_Up = Differentiation(theta_t_Up);
					Matrix<T> omega_t_Down = Differentiation(theta_t_Up);

					omega_phi_Up = omega_t_Up;
					omega_phi_Down = omega_t_Down;
					for (int i = 0; i < omega_phi_Up.Column(); i++)
					{
						omega_phi_Up[0][i] = phi_t_Up[1][i];
					}
					for (int i = 0; i < omega_phi_Down.Column(); i++)
					{
						omega_phi_Down[0][i] = phi_t_Down[1][i];
					}



				}

				static void MotorRotationalVelocityWithWalkingBeamAngle(T r, T R1, T R2, T H, T W, T L, std::function<T(T)> func, Matrix<T> &o_p_m, Matrix<T> &o_p_mDown, Multinomial<T> &omega_phi, Multinomial<T> &omega_phiDown)
				{

					SpatialBeamSolver<T>* spsolver = nullptr;
					SpatialBeamPostprocessor<T> sppost(spsolver);

					Multinomial<T> m = sppost.TaylorExpandWalkingBeamAngleToMotorRotationalVelocity(r, R1, R2, H, W, L);
					Multinomial<T> dm = m.Differentiate();
					Multinomial<T> mDown = sppost.TaylorExpandWalkingBeamAngleToMotorRotationalVelocity(r, R1, R2, H, W, L, false);
					Multinomial<T> dmDown = mDown.Differentiate();

					Matrix<T> m_Mul(1,m.n);
					for (int i = 0; i < m_Mul.Column(); i++)
					{
						m_Mul[0][i] = m.coefficients[i];
					}
					WriteMatrixToCSV("m_Mul.csv", m_Mul);

					std::function<T(T)> sPhi = func;

					Multinomial<T> sxnPhi = Taylor<T>(sPhi, 0, 3);
					Multinomial<T> sxndPhi = sxnPhi.Differentiate();

					Matrix<T> dPhi_Phi_points(2, 100);
					for (int i = 0; i < dPhi_Phi_points.Column(); i++)
					{
						dPhi_Phi_points[0][i] = sPhi(M_PI / 100 * i);
						dPhi_Phi_points[1][i] = NumericalDifferentiation<T>(M_PI / 100 * i, 1, sPhi, 0.0001);
					}
					Matrix<T> dPhi_Phi_Down_points(2, 100);
					for (int i = 0; i < dPhi_Phi_Down_points.Column(); i++)
					{
						dPhi_Phi_Down_points[0][i] = sPhi(M_PI + M_PI / 100 * i);
						dPhi_Phi_Down_points[1][i] = NumericalDifferentiation<T>(M_PI + M_PI / 100 * i, 1, sPhi, 0.0001);
					}
					WriteMatrixToCSV("dphi_phi_points.csv", dPhi_Phi_points);
					WriteMatrixToCSV("dphi_phi_down_points.csv", dPhi_Phi_Down_points);
					Multinomial<T> dPhi_Phi = LeastSquaresMethodToMultinomial(dPhi_Phi_points, 4);
					Multinomial<T> dPhi_Phi_Down = LeastSquaresMethodToMultinomial(dPhi_Phi_Down_points, 4);

					/*
					Matrix<T> dPhi_Phi_Up(2, 100);
					for (int i = 0; i < dPhi_Phi_Up.Column(); i++)
					{
					dPhi_Phi_Up[0][i] = 0.6979 + (2.3094 - 0.6979) / 100 * i;
					dPhi_Phi_Up[1][i] = dPhi_Phi_Up.Solve(o_p_mDown[0][i] - 1);
					}
					*/

					Matrix<T> dp_p_m(2, 100);
					for (int i = 0; i < dp_p_m.Column(); i++)
					{
						dp_p_m[0][i] = 0.88556 + (2.10221 - 0.88556) / 100 * i;
						dp_p_m[1][i] = dPhi_Phi.Solve(dp_p_m[0][i]);
					}

					Matrix<T> dp_p_mDown(2, 100);
					for (int i = 0; i < dp_p_mDown.Column(); i++)
					{
						dp_p_mDown[0][i] = 0.88556 + (2.10221 - 0.88556) / 100 * i;
						dp_p_mDown[1][i] = dPhi_Phi_Down.Solve(dp_p_mDown[0][i]);
					}
					WriteMatrixToCSV("dp_p.csv", dp_p_m);
					WriteMatrixToCSV("dp_pDown.csv", dp_p_mDown);
					//-std::asin((1-1.5037)/0.8057)
					//Multinomial<T> dPhi_Phi = Taylor<T>([&](T t) {return sxnPhi.Solve(t); }, [&](T t) {return sxndPhi.Solve(t); }, M_PI - std::acos((1 - 1.5037) / -0.8057),4);
					//Multinomial<T> dPhi_PhiDown = Taylor<T>([&](T t) {return sxnPhi.Solve(t); }, [&](T t) {return sxndPhi.Solve(t); }, std::acos((1 - 1.5037) / -0.8057),4);
					Multinomial<T> dm_ex = dm.X0Expand(-1);
					Multinomial<T> dmDown_ex = dmDown.X0Expand(-1);
					omega_phi = dm_ex.Multiply(dPhi_Phi);
					omega_phiDown = dmDown_ex.Multiply(dPhi_Phi_Down);

					o_p_m = Matrix<T>(2, 100);
					for (int i = 0; i < o_p_m.Column(); i++)
					{
						o_p_m[0][i] = 0.88556 + (2.10221 - 0.88556) / 100 * i;
						o_p_m[1][i] = omega_phi.Solve(o_p_m[0][i]);
					}

					o_p_mDown = Matrix<T>(2, 100);
					for (int i = 0; i < o_p_mDown.Column(); i++)
					{
						o_p_mDown[0][i] = 0.88556 + (2.10221 - 0.88556) / 100 * i;
						o_p_mDown[1][i] = omega_phiDown.Solve(o_p_mDown[0][i]);
					}

					std::function<T(T)> m_t = [&](T t) {
						return m.Solve(sPhi(t));
					};
					Matrix<T> m_t_m(2, 100);
					for (int i = 0; i < m_t_m.Column(); i++)
					{
						m_t_m[0][i] = 0 + M_PI / 100 * i;
						m_t_m[1][i] = m_t(m_t_m[0][i]);
					}

					Multinomial<T> m_t_T = Taylor<T>(m_t, 0, 4);
					Matrix<T> m_t_T_m(2, 100);
					for (int i = 0; i < m_t_T_m.Column(); i++)
					{
						m_t_T_m[0][i] = 0 + M_PI / 100 * i;
						m_t_T_m[1][i] = m_t_T.Solve(m_t_T_m[0][i] - 0);
					}

					Multinomial<T> mmm = Taylor<T>([&](T x) {return std::sin(x + 5); }, 1.5, 4);
					Matrix<T> mmm_m(2, 100);
					for (int i = 0; i < mmm_m.Column(); i++)
					{
						mmm_m[0][i] = 0 + M_PI / 100 * i;
						mmm_m[1][i] = mmm.Solve(mmm_m[0][i] - 1.5);
					}

					WriteMatrixToCSV("mmm_m.csv", mmm_m);
					WriteMatrixToCSV("m_t_T.csv", m_t_T_m);
					WriteMatrixToCSV("m_t.csv", m_t_m);
					WriteMatrixToCSV("o_p.csv", o_p_m);
					WriteMatrixToCSV("o_pDown.csv", o_p_mDown);

				}
			private:

			protected:

			};
		}
	}
}
