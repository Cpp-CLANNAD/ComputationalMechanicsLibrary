#pragma once
#include "../Core/Core.h"

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
			class PlaneRigidFrameSolver :public ISolver<T>
			{
			public:
				PlaneRigidFrameSolver(std::vector<IElement<T>*>& _element)
				{
					//assemb integral matrix
					this->element = _element;
					int length = this->Length();
					this->stiffness = Matrix<T>(length, length);
					this->force = Matrix<T>(length, 1);
					this->displacement = Matrix<T>(length, 1);
					this->unKnowForce.resize(length);
					this->unKnowDisplacement.resize(length);

					for (IElement<T>* e : this->element)
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

							this->force[(ni - 1) * 3 + i % 3][0] += e->Force()[i][0];
							this->unKnowForce[(ni - 1) * 3 + i % 3] = e->UnKnowForce()[i];
							this->displacement[(ni - 1) * 3 + i % 3][0] += e->Displacement()[i][0];
							this->unKnowDisplacement[(ni - 1) * 3 + i % 3] = e->UnKnowDisplacement()[i];


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

								this->stiffness[(ni - 1) * 3 + i % 3][(nj - 1) * 3 + j % 3] += e->Stiffness()[i][j];
							}
						}
				

					}
					
					//Cut Matrix: displacement known
					auto CutMatrix=[&](){
						//remove known displacements
						int newLength = 0;
						for (int f : this->unKnowDisplacement)
						{
							if (f == 0)
							{
								newLength++;
							}
						}

						Matrix<T> newStiffness(newLength, newLength);
						Matrix<T> newForce(newLength, 1);
						for (int i = 0, k = 0; i<newStiffness.Row(); k++, i++)
						{
							while (k < this->unKnowDisplacement.size() && this->unKnowDisplacement[k] == 1)
							{
								k++;
							}
							if (k >= this->unKnowDisplacement.size())
							{
								continue;
							}
							else
							{
								newForce[i][0] = this->force[k][0];
							}

							for (int j = 0, l = 0; j<newStiffness.Column(); j++, l++)
							{
								while (l < this->unKnowDisplacement.size() && this->unKnowDisplacement[l] == 1)
								{
									l++;
								}
								if (l >= this->unKnowDisplacement.size())
								{
									continue;
								}
								else
								{
									newStiffness[i][j] = this->stiffness[k][l];
								}

							}
						}

						Matrix<T> newDisplacement = (newStiffness^-1) * newForce;

						//full origin matrix
						for (int i = 0, j = 0; i<length; i++)
						{
							if (this->unKnowDisplacement[i] == 0)
							{
								this->displacement[i][0] = newDisplacement[j][0];
								j++;
							}
						}

						this->force = this->stiffness * this->displacement;
					};
					//Cut Zero Set One
					auto CutZero = [&](){

						int newLength = 0;
						for(int f : this->unKnowDisplacement)
						{
							if (f == 0)
							{
								newLength++;
							}
						}

						Matrix<T> newStiffness(this->stiffness);
						Matrix<T> newForce(this->force);

						for (int i = 0; i<this->unKnowDisplacement.size(); i++)
						{
							if (this->unKnowDisplacement[i] != 0)
							{
								newStiffness[i][i] = 1;
								newForce[i][0] = this->displacement[i][0];
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

						this->displacement = (newStiffness ^ -1) * newForce;

						this->force = this->stiffness * this->displacement;
					};
					
					//use it or not
                    ComputationalMechanicsLibrary::FiniteElementMethod::CutZeroAddOne(this->stiffness,this->displacement,this->force,this->unKnowDisplacement,this->unKnowForce,length);
				
				}

				std::vector<IElement<T>*>& Element()
				{
					return this->element;
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

				int Length()
				{
					//consecutive node number : 1,2,3,4,...
					int all = 0;
					for (IElement<T>* e : this->element)
					{
						for (int i : e->Node())
						{
							if (i > all)
							{
								all = i;
							}
						}
					}

					return all*3;
				}
			private:
				std::vector<IElement<T>*> element;
				Matrix<T> force;
				Matrix<T> displacement;
				Matrix<T> stiffness;
				std::vector<int> unKnowForce;
				std::vector<int> unKnowDisplacement;
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
				PlaneRigidFrameMaterial(T _elasticityModulus)
				{
					this->elasticityModulus = _elasticityModulus;
				}

				T ElasticityModulus()
				{
					return this->elasticityModulus;
				}
				T PoissonRatio()
				{
					return (T)0;
				}
				T Density()
				{
					return (T)0;
				}
			private:
				T elasticityModulus;
			protected:

			};
		}
	}
}

