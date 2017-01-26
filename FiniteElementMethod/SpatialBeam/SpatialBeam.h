#pragma once
#include "../Core/Core.h"

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
				SpatialBeamElement(int _i,int _j, IDynamicSection<T>& _section, double _length, double _angle, IDynamicMaterial<T>& _Material, std::vector<Matrix<T>>& _force, Matrix<T>& _displacement, Matrix<T> _velocity,Matrix<T> _acceleration)
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

					this->CalculateMass();
					this->CalculateStiffness();
				}

				std::vector<int>& Node()
				{
					return this->node;
				}

				std::vector<Matrix<T>>& Force()
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

				void DynamicReset()
				{
					this->CalculateStiffness();
				}
			private:
				std::vector<int> node;

				std::vector<Matrix<T>> force;

				Matrix<T> displacement;
				Matrix<T> velocity;
				Matrix<T> acceleration;

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
				T angle;

				void CalculateMass()
				{
					T l = this->length;
					T IxA = this->inertiaMomentX / this->sectionArea;
					this->mass = {{140, 0, 0, 0, 0, 0, 70, 0, 0, 0, 0, 0},
						{0, 156, 0, 0, 0, 22 * l, 0, 54, 0, 0, 0, -13 * l},
						{0, 0, 156, 0, -22 * l, 0, 0, 0, 54, 0, 13 * l, 0},
						{0, 0, 0, 40 * IxA, 0, 0, 0, 0, 0, 70 * IxA, 0, 0},
						{0, 0, -22 * l, 0, 4 * l*l, 0, 0, 0, -13 * l, 0, -3 * l*l, 0},
						{0, 22 * l, 0, 0, 0, 4 * l*l, 0, 13 * l, 0, 0, 0, -3 * l*l},
						{70, 0, 0, 0, 0, 0, 140, 0, 0, 0, 0, 0},
						{0, 54, 0, 0, 0, 13 * l, 0, 156, 0, 0, 0, -22 * l},
						{0, 0, 54, 0, -13 * l, 0, 0, 0, 156, 0, 22 * l, 0},
						{0, 0, 0, 70 * IxA, 0, 0, 0, 0, 0, 140 * IxA, 0, 0},
						{0, 0, 13 * l, 0, -3 * l*l, 0, 0, 0, 22 * l, 0, 4 * l*l, 0},
						{0, -13 * l, 0, 0, 0, -3 * l*l, 0, -22 * l, 0, 0, 0, 4 * l*l}};
					this->mass = (this->density*this->sectionArea*this->length / 420)* this->mass;
				}
				void CalculateStiffness()
				{
					T dix = this->displacement[0][0];
					T djx = this->displacement[6][0];
					T ��ix = this->displacement[3][0];
					T ��jx = this->displacement[9][0];
					T �� = this->poissonRatio;
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

					Matrix<T> Kl = {{E*A / l, 0, 0, 0, 0, 0, -E*A / l, 0, 0, 0, 0, 0},
						{0, 12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l, 0, -12 * E*Iz / l / l / l, 0, 0, 0, 6 * E*Iz / l / l},
						{0, 0, 12 * E*Iy / l / l / l, 1, -6 * E*Iy / l / l, 0, 0, 0, -12 * E*Iy / l / l / l, 0, -6 * E*Iy / l / l, 0},
						{0, 0, 0, G*Jx / l, 0, 0, 0, 0, 0, -G*Jx / l, 0, 0},
						{0, 0, -6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0},
						{0, 6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l, 0, 0, -6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l},
						{-E*A / l, 0, 0, 0, 0, 0, E*A / l, 0, 0, 0, 0, 0},
						{0, -12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l, 0, 12 * E*Iz / l / l / l, 0, 0, 0, -6 * E*Iz / l / l},
						{0, 0, -12 * E*Iy / l / l / l, 0, 6 * E*Iz / l / l, 0, 0, 12 * E*Iy / l / l / l, 0, 6 * E*Iy / l / l, 0},
						{0, 0, 0, -G*Jx / l, 0, 0, 0, 0, 0, G*Jx / l, 0, 0},
						{0, 0, -6 * E*Iy / l / l, 0, 2 * E*Iy / l, 0, 0, 0, 6 * E*Iy / l / l, 0, 4 * E*Iy / l, 0},
						{0, 6 * E*Iz / l / l, 0, 0, 0, 2 * E*Iz / l, 0, -6 * E*Iz / l / l, 0, 0, 0, 4 * E*Iz / l}};
					
					Matrix<T> KNB = {{3 / 2, 0, 0, 0, 0, 0, -3 / 2, 0, 0, 0, 0, 0},
						{0, 6 / 5, 0, 0, 0, l / 10, 0, -6 / 5, 0, 0, 0, l / 10},
						{0, 0, 6 / 5, 0, -l / 10, 0, 0, 0, -6 / 5, 0, -l / 10, 0},
						{0, 0, 0, Ix / A, 0, 0, 0, 0, 0, -Ix / A, 0, 0},
						{0, 0, -l / 10, 0, 2 * l*l / 15, 0, 0, 0, l / 10, 0, -l*l / 30, 0},
						{l / 10, 0, 0, 0, 0, 2 * l*l / 15, -l / 10, 0, 0, 0, 0, -l*l / 30},
						{-3 / 2, 0, 0, 0, 0, 0, 3 / 2, 0, 0, 0, 0, 0},
						{0, -6 / 5, 0, 0, 0, -l / 10, 0, 6 / 5, 0, 0, 0, -l / 10},
						{0, 0, -6 / 5, 0, l / 10, 0, 0, 0, 6 / 5, 0, l / 10, 0},
						{0, 0, 0, -Ix / A, 0, 0, 0, 0, 0, Ix / A, 0, 0},
						{0, 0, -l / 10, 0, -l*l / 30, 0, 0, 0, l / 10, 0, 2 * l*l / 15, 0},
						{l / 10, 0, 0, 0, 0, -l*l / 30, -l / 10, 0, 0, 0, 0, 2 * l*l / 15}};
					KNB = E*A*(djx - dix) / l / l*KNB;

					Matrix<T> KNT = {{0, 0, 0, �� / 2, 0, 0, 0, 0, 0, -�� / 2, 0, 0},
						{0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1, 0},
						{0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0, 1},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						{0, -1, 0, 0, 0, 0, 0, 1, 0, 0, 0, -l / 2},
						{0, 0, -1, 0, 0, 0, 0, 0, 1, 0, l / 2, 0},
						{0, 0, 0, -�� / 2, 0, 0, 0, 0, 0, �� / 2, 0, 0},
						{0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0},
						{0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -1},
						{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
						{0, 1, 0, 0, 0, l / 2, 0, -1, 0, 0, 0, 0},
						{0, 0, 1, 0, -�� / 2, 0, 0, 0, -1, 0, 0, 0}};
					KNT = (1 + 2 * ��)*G*Ix*(��jx - ��ix) / l / l * KNT;

					this->stiffness = Kl + KNB + KNT;
				}
			protected:

			};

			template<typename T>
			class SpatialBeamSolver :public IDynamicSolver<T>
			{
			public:
				SpatialBeamSolver(std::vector<IDynamicElement<T>*>& _element,T _timeInterval)
				{
					this->timeInterval = _timeInterval;
					//assemb integral matrix
					this->element = _element;
					int length = this->Length();

					this->mass=Matrix<T>(length, length);
					this->damp = Matrix<T>(length, length);
					this->stiffness = Matrix<T>(length, length);
					for (int i = 0; i < this->element[0]->Force().size(); i++)
					{
						this->force.push_back(Matrix<T>(length, 1));
						this->displacement.push_back(Matrix<T>(length, 1));
						this->velocity.push_back(Matrix<T>(length, 1));
						this->acceleration.push_back(Matrix<T>(length, 1));
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
							}

							this->displacement[0][(ni - 1) * 6 + i % 6][0] += e->Displacement()[i][0];
							this->velocity[0][(ni - 1) * 6 + i % 6][0] += e->Velocity()[i][0];
							this->acceleration[0][(ni - 1) * 6 + i % 6][0] += e->Acceleration()[i][0];


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
								//this->damp[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Damp()[i][j];
								this->stiffness[(ni - 1) * 6 + i % 6][(nj - 1) * 6 + j % 6] += e->Stiffness()[i][j];
							}
						}


					}

					for (int i = 0; i < this->force.size() - 1; i++)
					{
						this->Newmark(this->stiffness, this->mass, this->damp, this->displacement[i], this->velocity[i], this->acceleration[i], this->displacement[i + 1], this->velocity[i + 1], this->acceleration[i + 1], this->force[i + 1]);

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

								e->Displacement()[i][0] = this->displacement[0][(ni - 1) * 6 + i % 6][0];
							}
							//ReCalculate Stiffness nonliner
							e->DynamicReset();

							//full Stiffness
							this->stiffness = Matrix<T>(length, length);
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
								}
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
				std::vector<IDynamicElement<T>*> element;
				std::vector<Matrix<T>> force;

				std::vector<Matrix<T>> displacement;
				std::vector<Matrix<T>> velocity;
				std::vector<Matrix<T>> acceleration;

				Matrix<T> mass;
				Matrix<T> damp;
				Matrix<T> stiffness;

				T timeInterval;
				/// <summary>
				/// Nodes the iteration.
				/// kd=f
				/// </summary>
				/// <param name="K">The k.</param>
				/// <param name="d">The d.</param>
				/// <param name="F">The f.</param>
				void NodeIteration(Matrix<T>& K, Matrix<T>& d, Matrix<T>& F)
				{
					//Three diagonal matrix??
					d[1][0] = (F[0][0] - d[0][0] * K[0][0]) / K[0][1];
					for (int i = 2; i < d.Row(); i++)
					{
						d[i][0] = (F[i - 1][0] - d[i - 2][0] * K[i][i - 2] - d[i - 1][0] * K[i][i - 1]) / K[i - 1][i];
					}
				}
				void Newmark(Matrix<T>& K, Matrix<T>& M, Matrix<T>& C, Matrix<T>& d, Matrix<T>& v, Matrix<T>& a,Matrix<T>& dnext, Matrix<T>& vnext, Matrix<T>& anext , Matrix<T>& Fnext)
				{
					T DELTAt = this->timeInterval;
					T eta = 1.0 / 4;
					T delta = 1.0 / 2;
					T c[6] = { 1 / eta / DELTAt / DELTAt ,delta / eta / DELTAt / DELTAt ,1 / eta / DELTAt ,1 / 2 / eta - 1,delta / eta - 1, DELTAt*(delta / 2 / eta - 1) };
					Matrix<T> effectiveK = K + c[0] * M + c[1] * C;
					Matrix<T> effectiveFnext = Fnext + M*(c[0] * d + c[2] * v + c[3] * a) + C*(c[1] * d + c[4] * v + c[5] * a);
					
					//d(t+DELTAt),v(t+DELTAt),a(t+DELTAt)
					//this->NodeIteration(effectiveK, dnext, Fnext);
					dnext = effectiveK.Inverse()*Fnext;
					anext = c[0] * (dnext - d) - c[2] * v - c[3] * a;
					vnext = v + ((1 - delta)*a + delta*anext)*DELTAt;
				}

				void CalculateDamp()
				{
					T omega1=1, omega2=2;
					T zeta1=1, zeta2=2;

					T phi = 2 * (zeta1*omega2 - zeta2*omega1) / (omega2*omega2 - omega1*omega1);
					T xi = omega1*omega2*phi;

					this->damp = xi*this->mass + phi*this->stiffness;
				}
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