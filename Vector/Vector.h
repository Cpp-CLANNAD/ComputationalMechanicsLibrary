#pragma once
#include <vector>
#include <cmath>

namespace ComputationalMechanicsLibrary
{
	template<typename T>
	class Vector
	{
	public:

		enum VectorType
		{
			Row = 0,
			Column = 1
		};

		Vector(const VectorType _type = VectorType::Row)
		{
			this->array = nullptr;
			this->length = 0;
			this->type = _type;
		}

		Vector(const size_t _length,const VectorType _type=VectorType::Row)
		{
			this->array = new T[_length];
			this->length = _length;
			this->type = _type;
			for (int i = 0; i < this->length; i++)
			{
				this->array[i] = T(0);
			}
		}

		Vector(const Vector<T>& _vector)
		{
			this->FullArray(_vector);
		}

		Vector(const T _array[], const size_t _length,const VectorType _type = VectorType::Row)
		{
			this->array = new T[_length];
			this->length = _length;
			this->type = _type;
			for (int i = 0; i < this->length; i++)
			{
				this->array[i] = _array[i];
			}
		}

		Vector(const std::vector<T> _array, const VectorType _type = VectorType::Row)
		{
			this->array = new T[_array.size()];
			this->length = _array.size();
			this->type = _type;
			for (int i = 0; i < this->length; i++)
			{
				this->array[i] = _array[i];
			}
		}

		Vector(const std::initializer_list<T> _array, const VectorType _type = VectorType::Row)
		{
			this->array = new T[_array.size()];
			this->length = _array.size();
			this->type = _type;
			for (int i = 0; i < this->length; i++)
			{
				this->array[i] = *(_array.begin() + i);
			}
		}

		T* Array() const
		{
			return this->array;
		}

		VectorType Type() const
		{
			return this->type;
		}

		size_t Length() const
		{
			return (this->length);
		}

		~Vector()
		{
			delete[] this->array;
			this->array = nullptr;
			this->length = 0;
		}


		Vector<T> CreateRowTransportVector()
		{
			Vector<T> result((Vector<T>&)(*this));
			result.type = VectorType::Row;
			return result;
		}

		Vector<T> CreateColumnTransportVector()
		{
			Vector<T> result((Vector<T>&)(*this));
			result.type = VectorType::Column;
			return result;
		}

		Vector<T> Transport() const
		{
			Vector<T> result=*this;
			if (result.type == VectorType::Row)
			{
				result.type = VectorType::Column;
			}
			else if (result.type == VectorType::Column)
			{
				result.type = VectorType::Row;
			}
			return result;
		}

		static Vector<T> Transport(const Vector<T>& _vector)
		{
			return _vector.Transport();
		}

		Vector<T> Add(const Vector<T>& _vector1)
		{
			if (_vector1.type != this->type)
			{
				throw std::exception("type not equal");
			}
			if (_vector1.length != this->length)
			{
				throw std::exception("dimension not equal");
			}

			Vector<T> result(_vector1);
			for (size_t i = 0; i < result.length; i++)
			{
				result.array[i] += this->array[i];
			}

			return result;
		}

		Vector<T> Subtract(const Vector<T>& _vector1)
		{
			if (_vector1.type != this->type)
			{
				throw std::exception("type not equal");
			}
			if (_vector1.length != this->length)
			{
				throw std::exception("dimension not equal");
			}

			Vector<T> result(_vector1);
			for (size_t i = 0; i < result.length; i++)
			{
				result.array[i] -= this->array[i];
			}

			return result;
		}

		Vector<T> CrossProduct(const Vector<T>& _vector1)
		{
			/*
			if (_vector1.type != this->type)
			{
			throw std::exception("type not equal");
			}
			*/
			if (_vector1.length != 3 || this->length != 3)
			{
				throw std::exception("dimension not 3");
			}
			if (_vector1.length != this->length)
			{
				throw std::exception("dimension not equal");
			}

			Vector<T> result = { this->array[1] * _vector1.array[2] - this->array[2] * _vector1.array[1] ,this->array[2] * _vector1.array[0] - this->array[0] * _vector1.array[2],this->array[0] * _vector1.array[1] - this->array[1] * _vector1.array[0] };

			return result;

		}
		static Vector<T> CrossProduct(const Vector<T>& _vector1, const Vector<T>& _vector2)
		{
			return _vector1.CrossProduct(_vector2);
		}
		T ScalarProduct(const Vector<T>& _vector1)
		{
			/*
			if (_vector1.type != this->type)
			{
				throw std::exception("type not equal");
			}	
			*/

			if (_vector1.length != this->length)
			{
				throw std::exception("dimension not equal");
			}

			T result = T(0);

			for (size_t i = 0; i < this->length; i++)
			{
				result += this->array[i] * _vector1.array[i];
			}


			return result;
		}

		T Module()
		{
			T result = 0;
			for (size_t i = 0; i < this->length; i++)
			{
				result += std::pow(this->array[i]);
			}
			result = std::sqrt(result);

			return result;
		}

		T IncludedAngle(const Vector<T>& _vector1)
		{
			return std::acos(this->ScalarProduct(_vector1) / this->Module() / _vector1.Module());
		}
		static T IncludedAngle(const Vector<T>& _vector1,const Vector<T>& _vector2)
		{
			return _vector1.IncludedAngle(_vector2);
		}

		std::vector<T> Angle()
		{
			std::vector<T> result;
			result.reserve(this->length);

			for (size_t i = 0; i < this->length; i++)
			{
				Vector<T> il(this->length,this->type);
				il.array[i] = 1;

				result[i] = this->IncludedAngle(il);
			}

			return result;
		}

		Vector<T>& operator =(const Vector<T>& _vector) const
		{
			this->FullArray(_vector);
			return *this;
		}
		
		bool operator == (const Vector<T>& _vector) const
		{
			if (this->type != _vector.type)
			{
				return false;
			}
			if (this->length != _vector.length)
			{
				return false;
			}
			for (size_t i = 0; i < this->length; i++)
			{
				if (this->array[i] != _vector[i])
				{
					return false;
				}
			}

			return true;
		}

		friend Vector<T> operator + (Vector &a, const Vector<T> &b)
		{
			return a.Add(b);
		}

		friend Vector<T> operator - (Vector &a, const Vector<T> &b)
		{
			return a.Subtract(b);
		}

		friend Vector<T> operator * (Vector &a, const Vector<T> &b)
		{
			return a.ScalarProduct(b);
		}

	protected:

	private:
		T* array;
		VectorType type;
		size_t length;


		void FullArray(const ComputationalMechanicsLibrary::Vector<T> & _vector)
		{
			this->array = new T[_vector.length];
			this->length = _vector.length;
			this->type = _vector.type;
			for (int i = 0; i < this->length; i++)
			{
				this->array[i] = _vector.array[i];
			}
		}
	};

	template<typename T>
	Vector<T> Transport(const Vector<T>& _vector)
	{
		return _vector.Transport();
	}

	template<typename T>
	Vector<T> CrossProduct(const Vector<T>& _vector1, const Vector<T>& _vector2)
	{
		return _vector1.CrossProduct(_vector2);
	}

	template<typename T>
	T IncludedAngle(const Vector<T>& _vector1, const Vector<T>& _vector2)
	{
		return _vector1.IncludedAngle(_vector2);
	}

}