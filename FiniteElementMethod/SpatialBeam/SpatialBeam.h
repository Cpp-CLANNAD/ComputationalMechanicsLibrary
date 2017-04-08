#pragma once
#include "../Core/Core.h"

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
				SpatialBeamSolver(std::vector<SpatialBeamElement<T>*>& _element, T _timeInterval, DampType _type = DampType::Rayleigh,
					T _omega1 = 0.2507,T _omega2 = 0.4413,T _zeta1 = 0.2,T _zeta2 = 0.2)
				{
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
					for (int t = 0; t < this->force.size() - 1; t++)
					{
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

					this->NodeIteration(effectiveK, dnext, effectiveFnext);

					Matrix<T> effectiveKv = 1.0 / delta / DELTAt*M + C;
					Matrix<T> effectiveFnextv = Fnext + M*(1.0 / delta/DELTAt*v + (1.0 / delta + 1)*a) - K*dnext;

					WriteMatrixToCSV("effKv.csv", effectiveKv);
					WriteMatrixToCSV("effFv.csv", effectiveFnextv);

					this->NodeIteration(effectiveKv, vnext, effectiveFnextv);

					//WriteMatrixToCSV("effVnext.csv", vnext);

					Matrix<T> effectiveKa = M;
					Matrix<T> effectiveFnexta = F - C*vnext - K*dnext;

					WriteMatrixToCSV("effKa.csv", effectiveKa);
					WriteMatrixToCSV("effFa.csv", effectiveFnexta);

					this->NodeIteration(effectiveKa, anext, effectiveFnexta);

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
					//////////////////////////anext = c[0] * (dnext - d) - c[2] * v - c[3] * a;
					//temp=v + ((1 - delta)*a + delta*anext)*DELTAt;
					//CutZeroAddOne(EM,vnext,temp,unKnowVnext,std::vector<int>(this->Length(),1),this->Length());
					/////////////////////////vnext = v + c[6] * a + c[7] * anext;
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

		}
	}
}
