#pragma once
#include <cstdlib>
#include <string>
#include <exception>
#include <initializer_list>
#include <vector>

#include <iostream>
#include <fstream>

#include "../Vector/Vector.h"

namespace ComputationalMechanicsLibrary
{
	/// <summary>
	/// the class of Matrix
	/// </summary>
	/// <typename name="T">the type of matrix's number</typename>
	template<typename T>
	class Matrix
	{
	public:

		/// <summary>
		/// for Inverse() parameter
		/// </summary>
		enum InverseMethod
		{			
			/// <summary>
			/// The gauss jordan
			/// </summary>
			GaussJordan = 1
		};

		/// <summary>
		/// Initializes a new instance of a ?null? Matrix(0x0) class.
		/// </summary>
		Matrix()
		{
			this->array = nullptr;
			this->row = 0;
			this->column = 0;
		}
		/// <summary>
		/// Initializes a new instance of the Matrix(1x1) class by a number.
		/// </summary>
		/// <param name="n">The number.</param>
		Matrix(T n)
		{
			this->array = new T*[1];
			this->array[0] = new T[1];
			this->array[0][0] = n;
			this->row = 1;
			this->column = 1;
		}
		/// <summary>
		/// Initializes a new instance of the Matrix class by a existed matrix.
		/// just like copy
		/// </summary>
		/// <param name="a">existed matrix</param>
		Matrix(const Matrix<T>& a)
		{
			this->array = nullptr;
			this->FullArray(a);
		}
		/// <summary>
		/// Initializes a new instance of the matrix class by set row and column
		/// </summary>
		/// <param name="_row">The row.</param>
		/// <param name="_column">The column.</param>
		Matrix(size_t _row, size_t _column)
		{
            this->array = new T*[_row];
			this->array[0] = new T[_row*_column]{ 0 };
			for (size_t i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
			}
			this->row = _row;
			this->column = _column;
		}
		/// <summary>
		/// Initializes a new instance of the Matrix class by array
		/// </summary>
		/// <param name="arr">The arr.</param>
		/// <param name="_row">The array's row.</param>
		/// <param name="_column">The array's column.</param>
		Matrix(T** arr, size_t _row, size_t _column)
		{
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (size_t j = 0; j < _column; j++)
			{
				this->array[0][j] = arr[0][j];
			}
			for (size_t i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (size_t j = 0; j < _column; j++)
				{
					this->array[i][j] = arr[i][j];
				}
			}
			this->row = _row;
			this->column = _column;
		}
		/// <summary>
		/// Initializes a new instance of the <see cref="Matrix"/> class.
		/// </summary>
		/// <param name="arr">The arr.</param>
		Matrix(std::vector<std::vector<T>> arr)
		{
			size_t _row = arr.size();
			size_t _column = arr[0].size();
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (size_t j = 0; j < _column; j++)
			{
				this->array[0][j] = arr[0][j];
			}
			for (size_t i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (size_t j = 0; j < _column; j++)
				{
					this->array[i][j] = arr[i][j];
				}
			}
			this->row = _row;
			this->column = _column;
		}
		/// <summary>
		/// Initializes a new instance of the <see cref="Matrix"/> class.
		/// </summary>
		/// <param name="arr">The arr.</param>
		Matrix(std::initializer_list<std::initializer_list<T>> arr)
		{
			size_t _row = arr.size();
			size_t _column = (arr.begin())->size();
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (size_t j = 0; j < _column; j++)
			{
				this->array[0][j] = *((arr.begin())->begin()+j);
			}
			for (size_t i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (size_t j = 0; j < _column; j++)
				{
					this->array[i][j] = *((arr.begin()+i)->begin() + j);
				}
			}
			this->row = _row;
			this->column = _column;
		}

		Matrix<T>(const Vector<T>& _vector)
		{
			size_t _row = (_vector.Length() - 1)*_vector.Type() + 1;
			size_t _column = _vector.Length()/_row;
			this->array = new T*[_row];
			this->array[0] = new T[_vector.Length()]{ 0 };
			for (size_t i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
			}

			size_t t = 0;
			for (size_t i = 0; i < _row; i++)
			{
				for (size_t j = 0; j < _column; j++)
				{
					this->array[i][j] = _vector.Array()[t];
					t++;
				}
			}

			this->row = _row;
			this->column = _column;
		}

		/// <summary>
		/// Finalizes an instance of the Matrix class.
		/// </summary>
		~Matrix()
		{
			this->ClearArray();
		}

        void SetUnit(size_t _row)
        {
            this->destory();
            this->array = new T*[_row];
            this->array[0] = new T[_row*_row]{0};
            for (size_t j = 0; j < _row; j++)
            {
                this->array[0][j] = ((T*)this->array)[0 * _row + j];
            }
            for(size_t i=0;i<_row;i++)
            {
                this->array[i][i]=1;
            }
            this->row = _row;
            this->column = _row;
        }
		
		//full this matrix by another matrix, just change values in array, but not change address.
		void FullMatrix(Matrix<T> &m)
		{
			if (m.Row() != this->row || m.Column() != this->column)
			{
				throw std::exception("not align row or column");
			}

			for (size_t i = 0; i < this->row; i++)
			{
				for (size_t j = 0; j < this->column; j++)
				{
					this->array[i][j] = m[i][j];
				}
			}
		}

		/// <summary>
		/// Check whether this matrix is Square matrix
		/// </summary>
		/// <returns>Square is true, or else is false</returns>
		bool CheckSquare()
		{
			return this->row == this->column ? true : false;
		}
		/// <summary>
		/// Get the row.
		/// </summary>
		/// <returns>row</returns>
		size_t Row() const
		{
			return this->row;
		}
		/// <summary>
		/// Get the Column.
		/// </summary>
		/// <returns>column</returns>
		size_t Column() const
		{
			return this->column;
		}
		/// <summary>
		/// Get the Array which is used to store the Matrix
		/// </summary>
		/// <returns>A point, point to the array</returns>
		T** Array() const
		{
			return this->array;
		}
		
		/// <summary>
		/// Ats the specified array[i][j]
		/// </summary>
		/// <param name="i">The i.</param>
		/// <param name="j">The j.</param>
		/// <returns></returns>
		T at(size_t i, size_t j)
		{
			if (i < 0 || i >= this->row || j < 0 || j >= this->column)
			{
				throw std::exception("out matrix length");
			}

			return this->array[i][j];
		}

		/// <summary>
		/// Determinant
		/// </summary>
		/// <param name="method">The method.</param>
		/// <returns>this matrix's determinat</returns>
		T Determinant(InverseMethod method = InverseMethod::GaussJordan)
		{
			T result;
			try
			{
				switch (method)
				{
				case InverseMethod::GaussJordan:
					result = this->DeterminantGaussJordan();
					break;
				default:
					result = this->DeterminantGaussJordan();
					break;
				}

				return result;
			}
			catch (...)
			{
				throw;
			}
		}
		/// <summary>
		/// Inverse
		/// </summary>
		/// <param name="method">The method.</param>
		/// <returns>this matrix's inverse matrix</returns>
		Matrix<T> Inverse(InverseMethod method = InverseMethod::GaussJordan)
		{
			Matrix<T> result;
			try
			{
				switch (method)
				{
				case InverseMethod::GaussJordan:
					result = this->InverseGaussJordan();
					break;
				default:
					result = this->InverseGaussJordan();
					break;
				}

				return result;
			}
            catch (std::exception &e)
			{
				throw;
			}
		}
		/// <summary>
		/// Transport
		/// </summary>
		/// <returns>this matrix's transport matrix</returns>
		Matrix<T> Transport()
		{
			Matrix<T> result(this->column, this->row);
			for (size_t i = 0; i < result.row; i++)
			{
				for (size_t j = 0; j < result.column; j++)
				{
					result.array[i][j] = this->array[j][i];
				}
			}
			return result;
		}
		/// <summary>
		/// Adjoint
		/// </summary>
		/// <returns>this matrix's adjoint matrix</returns>
		Matrix<T> AdjointMatrix()
		{
			//阶数
			size_t n = this->CheckSquare() ? this->row : throw std::exception("not a square");
			//顺便把异常抛了出去
			//制作一个伴随矩阵大小的矩阵
			Matrix<T> result(n, n);

			//存储代数余子式的矩阵（行、列数都比原矩阵少1）
			Matrix<T> temp(n - 1, n - 1);
			//生成伴随矩阵
			for (size_t i = 0; i < n; i++)
			{
				for (size_t j = 0; j < n; j++)
				{
					//生成代数余子式
					for (size_t x = 0; x < n-1; x++)
					{
						for (size_t y = 0; y < n-1; y++)
						{
							temp.array[x][y] = this->array[x < i ? x : x + 1][y < j ? y : y + 1];
						}
					}
					result.array[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * temp.Determinant();
				}
			}

			return result;
		}

		/// <summary>
		/// Add with other matrix
		/// </summary>
		/// <param name="b">addend matrix</param>
		/// <returns>sum matrix</returns>
		Matrix<T> Add(const Matrix<T> &b)
		{
			if (this->CheckDimension(b)==false)
			{
				throw std::exception("dimensions not equal");
			}
			Matrix<T> result(*this);
			for (size_t i = 0; i < result.row; i++)
			{
				for (size_t j = 0; j < result.column; j++)
				{
					result.array[i][j] += b.array[i][j];
				}
			}

			return result;
		}
		/// <summary>
		/// Subtract with other matrix
		/// </summary>
		/// <param name="b">subtrahend matrix</param>
		/// <returns>difference matrix</returns>
		Matrix<T> Subtract(const Matrix<T> &b)
		{
			if (this->CheckDimension(b)==false)
			{
				throw std::exception("dimensions not equal");
			}
			Matrix<T> result(*this);
			for (size_t i = 0; i < result.row; i++)
			{
				for (size_t j = 0; j < result.column; j++)
				{
					result.array[i][j] -= b.array[i][j];
				}
			}

			return result;
		}
		/// <summary>
		/// Multiplie with other matrix
		/// </summary>
		/// <param name="b">multiplier matrix</param>
		/// <returns>product matrix</returns>
		virtual Matrix<T> Multiply(const Matrix<T> &b)
		{
			if (this->column != b.row)
			{
				throw std::exception("dimension not alignment");
			}
			Matrix result(this->row, b.column);
			for (size_t i = 0; i < result.row; i++)
			{
				for (size_t j = 0; j < result.column; j++)
				{
					for (size_t k = 0; k < this->column; k++)
					{
						result.array[i][j] += this->array[i][k] * b.array[k][j];
					}
				}
			}
			return result;
		}

		Matrix<T> GeneralizedInverse()
		{
			if (this->row == this->column)
			{
				return this->Inverse();
			}
			else if (this->row > this->column)
			{
				Matrix<T> A(this->column, this->column);
				Matrix<T> B(this->row - this->column, this->column);
				for (size_t i = 0; i < this->row; i++)
				{
					for (size_t j = 0; j < this->column; j++)
					{
						if (i < this->column)
						{
							A[i][j] = this->array[i][j];
						}
						else
						{
							B[i-this->column][j] = this->array[i][j];
						}
					}
				}

				Matrix<T> E(this->column,this->column);
				for (size_t i = 0; i < E.Row(); i++)
				{
					E[i][i] = 1;
				}
				WriteMatrixToCSV("A.csv", A);
				Matrix<T> BAi = B*(A^-1);

				Matrix<T> D = A;
				Matrix<T> C(this->row, this->column);
				for (size_t i = 0; i < this->row; i++)
				{
					for (size_t j = 0; j < this->column; j++)
					{
						if (i < this->column)
						{
							C[i][j] = E[i][j];
						}
						else
						{
							C[i][j] = BAi[i - this->column][j];
						}
					}
				}

				Matrix<T> result = D.Transport()*((D*D.Transport()) ^ -1)*((C.Transport()*C) ^ -1)*C.Transport();
				return result;

			}
			else
			{
				Matrix<T> A(this->row, this->row);
				Matrix<T> B(this->row, this->column - this->row);
				for (size_t i = 0; i < this->row; i++)
				{
					for (size_t j = 0; j < this->column; j++)
					{
						if (j < this->row)
						{
							A[i][j] = this->array[i][j];
						}
						else
						{
							B[i][j- this->row] = this->array[i][j];
						}
					}
				}

				Matrix<T> E(this->column, this->column);
				for (size_t i = 0; i < E.Row(); i++)
				{
					E[i][i] = 1;
				}
				Matrix<T> AiB = (A^-1)*B;

				Matrix<T> C = A;
				Matrix<T> D(this->row, this->column);
				for (size_t i = 0; i < this->row; i++)
				{
					for (size_t j = 0; j < this->column; j++)
					{
						if (j < this->row)
						{
							D[i][j] = E[i][j];
						}
						else
						{
							D[i][j] = AiB[i][j - this->row];
						}
					}
				}
				Matrix<T> result = D.Transport()*((D*D.Transport()) ^ -1)*((C.Transport()*C) ^ -1)*C.Transport();
				return result;
			}
		}

		/// <summary>
		/// Check the Dimensions whether is different from other matrix
		/// </summary>
		/// <param name="b">comparison matrix</param>
		/// <returns>the same is true, or else is false</returns>
		bool CheckDimension(const Matrix<T> &b)
		{
			return (this->row == b.row) && (this->column == b.column);
		}

		//operator

		/// <summary>
		/// Operator=s the specified b.
		/// </summary>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		Matrix<T>& operator = (const Matrix<T> &b) const
		{
			this->FullArray(b);
			return *this;
		}
		/// <summary>
		/// Operator==s the specified b.
		/// </summary>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		bool operator == (const Matrix<T> &b) const
		{
			if (this->row != b.row || this->column != b.column)
			{
				if (this->array == nullptr&&b.array == nullptr)
				{
					return true;
				}
				else if (this->array != nullptr&&b.array != nullptr)
				{
					for (size_t i = 0; i < this->row; i++)
					{
						for (size_t j = 0; j < this->column; j++)
						{
							if (this->array[i][j] != b.array[i][j])
							{
								return false;
							}
						}
					}
					return true;
				}
			}
			return false;
		}
		/// <summary>
		/// Operator+s the specified a.
		/// </summary>
		/// <param name="a">a.</param>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		friend Matrix<T> operator + (Matrix<T> &a, const Matrix<T> &b)
		{
			return a.Add(b);
		}
		/// <summary>
		/// Operator-s the specified a.
		/// </summary>
		/// <param name="a">a.</param>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		friend Matrix<T> operator - (Matrix<T> &a, const Matrix<T> &b)
		{
			return a.Subtract(b);
		}
		/// <summary>
		/// Operator*s the specified a.
		/// </summary>
		/// <param name="a">a.</param>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		friend Matrix<T> operator * (Matrix<T> &a, const Matrix<T> &b)
		{
			return a.Multiply(b);
		}
		/// <summary>
		/// Operator^s the specified a.
		/// </summary>
		/// <param name="a">a.</param>
		/// <param name="n">The n.</param>
		/// <returns></returns>
		friend Matrix<T> operator ^(Matrix<T> &a, int n)
		{
			Matrix<T> result;
			if (n > 0)
			{
				result = a;
				for (size_t i = 1; i < n; i++)
				{
					result = result.Multiply(a);
				}

			}
			else if (n == 0)
			{
				result =Matrix<T>(a.Row(), a.Column());
				for (size_t i = 0; i < result.Row(); i++)
				{
					result.Array()[i][i] = 1;
				}
			}
			else
			{
				Matrix<T> ainv = a.Inverse();
				result = Matrix<T>(ainv);
				for (size_t i = -1; i > n; i--)
				{
					result = result.Multiply(ainv);
				}
			}


			return result;
		}
		
		/// <summary>
		/// Operator*s the specified n.(number*matrix)
		/// </summary>
		/// <param name="n">The n.</param>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		friend Matrix<T> operator * (T n, Matrix<T> &b)
		{
			return b*n;
		}
		/// <summary>
		/// Operator*s the specified b.(matrix*number)
		/// </summary>
		/// <param name="b">The b.</param>
		/// <param name="n">The n.</param>
		/// <returns></returns>
		friend Matrix<T> operator * (Matrix<T> &b, T n)
		{
			Matrix<T> result(b);
			for (size_t i = 0; i < b.Row(); i++)
			{
				for (size_t j = 0; j < b.Column(); j++)
				{
					result.Array()[i][j] *= n;
				}
			}
			return result;
		}

		friend Matrix<T> operator *(Vector<T> &a, const Matrix<T> &b)
		{
			return Matrix<T>(a)*b;
		}

		/// <summary>
		/// matrix[i][j]===>>>matrix.array[i][j]
		/// </summary>
		/// <param name="i">The i.</param>
		/// <returns></returns>
		T* operator [] (size_t i)
		{
			return this->array[i];
		}

		/// <summary>
		/// destories this instance.
		/// for delete it
		/// matrix can't build in stack
		/// </summary>
		void destory()
		{
			delete this;
		}
		/// <summary>
		/// Destories the specified matrix
		/// for delete it
		/// matrix can't build in stack
		/// </summary>
		/// <param name="a">the matrix to destory</param>
		static void destory(Matrix<T>* a)
		{
			delete a;
		}
	private:
		/// <summary>
		/// The row
		/// </summary>
		size_t row;
		/// <summary>
		/// The column
		/// </summary>
		size_t column;
		/// <summary>
		/// point to the store array
		/// </summary>
		T** array;

		/// <summary>
		/// Clears the array.
		/// </summary>
		void ClearArray()
		{
			if (this->array != nullptr)
			{
				delete[] this->array[0];
				this->array[0] = nullptr;
				delete[] this->array;
				this->array = nullptr;
			}
			this->row = 0;
			this->column = 0;
		}
		/// <summary>
		/// Fulls the array by a existed matrix
		/// </summary>
		/// <param name="a">existed matrix</param>
        void FullArray(const Matrix<T> &a)
		{
			this->ClearArray();
            this->array = new T*[a.row];
			this->array[0] = new T[(a.row)*(a.column)];
			for (size_t j = 0; j < a.column; j++)
			{
				this->array[0][j] = a.array[0][j];
			}
			for (size_t i = 1; i < a.row; i++)
			{
				this->array[i] = (this->array[i - 1]) + (a.column);
				for (size_t j = 0; j < a.column; j++)
				{
					this->array[i][j] = a.array[i][j];
				}
			}
			this->row = a.row;
			this->column = a.column;
		}
		/// <summary>
		/// Swaps the Row of this matrix
		/// </summary>
		/// <param name="i">one row</param>
		/// <param name="j">another row</param>
		void SwapRowSelf(size_t i, size_t j)
		{
			T temp;
			for (size_t k = 0; k < this->column; k++)
			{
				temp = this->array[i][k];
				this->array[i][k] = this->array[j][k];
				this->array[j][k] = temp;
			}
		}
		/// <summary>
		/// Swaps the Column of this matrix
		/// </summary>
		/// <param name="i">one column</param>
		/// <param name="j">another column</param>
		void SwapColumnSelf(size_t i, size_t j)
		{
			T temp;
			for (size_t k = 0; k < this->row; k++)
			{
				temp = this->array[k][i];
				this->array[k][i] = this->array[k][j];
				this->array[k][j] = temp;
			}
		}
		/// <summary>
		/// Transport this matrix
		/// </summary>
		void TransportSelf()
		{
			Matrix* tempMatrix = this->Transport();
			this->FullArray(tempMatrix);
			delete tempMatrix;
		}

		/// <summary>
		/// Inverse by the Gauss Jordan method
		/// </summary>
		/// <returns>this matrix's inverse matrix</returns>
		Matrix<T> InverseGaussJordan()
		{
			//阶数
			size_t n = this->CheckSquare() ? this->row : throw std::exception("not a square");
			//顺便把异常抛了出去

			Matrix<T> result(*this);
			size_t* jr = new size_t[this->row];
			size_t* jc = new size_t[this->column];

			double max;

			for (size_t k = 0; k < result.row; k++)
			{
				//全选主元
				max = 0;
				jr[k] = k;
				jc[k] = k;

				for (size_t i = k; i < result.row; i++)
				{
					for (size_t j = k; j < result.column; j++)
					{
						if (max < std::abs(result.array[i][j]))
						{
							max = std::abs(result.array[i][j]);
							jr[k] = i;
							jc[k] = j;
						}
					}
				}

				if (max < 1.0e-24)
				{   //! 无逆矩阵
					delete[] jr;
					delete[] jc;
					throw std::exception("determinat is 0");
				}

				//交换
				if (jr[k] != k)
				{
					result.SwapRowSelf(k, jr[k]);
				}
				if (jc[k] != k)
				{
					result.SwapColumnSelf(k, jc[k]);
				}

				result.array[k][k] = 1.0 / result.array[k][k];

				for (size_t j = 0; j < result.column; j++)
				{
					if (j != k)
					{
						result.array[k][j] *= result.array[k][k];
					}
				}
				for (size_t i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						for (size_t j = 0; j < result.column; j++)
						{
							if (j != k)
							{
								result.array[i][j] -= result.array[i][k] * result.array[k][j];
							}
						}
					}
				}
				for (size_t i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						result.array[i][k] *= -result.array[k][k];
					}
				}

			}

			//恢复
			for (size_t k = result.row - 1; k >= 0; k--)
			{
				if (jc[k] != k)
				{
					result.SwapRowSelf(k, jc[k]);
				}
				if (jr[k] != k)
				{
					result.SwapColumnSelf(k, jr[k]);
				}
			}

			delete[] jr;
			delete[] jc;

			return result;
		}
		/// <summary>
		/// Determinant by the Gauss Jordan method
		/// </summary>
		/// <returns>
		/// this matrix's determinat
		/// </returns>
		T DeterminantGaussJordan()
		{
			//阶数
			size_t n = this->CheckSquare() ? this->row : throw std::exception("not a square");			
			//顺便把异常抛了出去
			T resultD = (T)1;


			Matrix<T> result(*this);
			size_t* jr = new size_t[this->row];
			size_t* jc = new size_t[this->column];

			double max;

			for (size_t k = 0; k < result.row; k++)
			{
				//全选主元
				max = 0;
				jr[k] = k;
				jc[k] = k;

				for (size_t i = k; i < result.row; i++)
				{
					for (size_t j = k; j < result.column; j++)
					{
						if (max < std::abs(result.array[i][j]))
						{
							max = std::abs(result.array[i][j]);
							jr[k] = i;
							jc[k] = j;
						}
					}
				}

				if (max < 1.0e-12)
				{   //! 无逆矩阵
					delete[] jr;
					delete[] jc;
					return (T)0;
				}

				//交换
				if (jr[k] != k)
				{
					resultD *= (T)-1;
					result.SwapRowSelf(k, jr[k]);
				}
				if (jc[k] != k)
				{
					resultD *= (T)-1;
					result.SwapColumnSelf(k, jc[k]);
				}

				resultD *= result.array[k][k];
				result.array[k][k] = 1.0 / result.array[k][k];

				for (size_t j = 0; j < result.column; j++)
				{
					if (j != k)
					{
						result.array[k][j] *= result.array[k][k];
					}
				}
				for (size_t i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						for (size_t j = 0; j < result.column; j++)
						{
							if (j != k)
							{
								result.array[i][j] -= result.array[i][k] * result.array[k][j];
							}
						}
					}
				}
				for (size_t i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						result.array[i][k] *= -result.array[k][k];
					}
				}

			}

			//恢复
			for (size_t k = result.row - 1; k >= 0; k--)
			{
				if (jc[k] != k)
				{
					result.SwapRowSelf(k, jc[k]);
				}
				if (jr[k] != k)
				{
					result.SwapColumnSelf(k, jr[k]);
				}
			}

			delete[] jr;
			delete[] jc;
			return resultD;

		}


		
	protected:

	};

	//Now the CSV file is stdand , the last in line(row) crlf without comma
	//Please don't use quote

	template<typename T>
	void WriteMatrixToCSV(std::string filename, ComputationalMechanicsLibrary::Matrix<T> &matrix)
	{
		std::ofstream outfile(filename);

		for (size_t i = 0; i < matrix.Row(); i++)
		{
			for (size_t j = 0; j < matrix.Column()-1; j++)
			{
				outfile << matrix[i][j] << ",";
			}
			outfile << matrix[i][matrix.Column() - 1];
			outfile << "\n";
		}

		outfile.close();
	}

	template<typename T>
	void ReadMatrixFromCSV(std::string filename, ComputationalMechanicsLibrary::Matrix<T> &matrix)
	{
		std::ifstream infile(filename);

		
		std::ifstream inFileRow(filename);
		std::string  s;
		size_t row = 0;
		size_t column = 0;
		bool ifColumn = false;
		while (std::getline(inFileRow, s))
		{
			row++;
			if (ifColumn == false)
			{
				column = 0;
				for (size_t i = 0; i < s.length(); i++)
				{
					if (s[i] == ',')
					{
						column++;
					}
				}
				column++;

				ifColumn = true;
			}

		}


		matrix = Matrix<T>(row, column);

		char d;
		for (size_t i = 0; i < matrix.Row(); i++)
		{
			for (size_t j = 0; j < matrix.Column()-1; j++)
			{
				infile >> matrix[i][j] >> d;
			}
			infile >> matrix[i][matrix.Column() - 1];
			//infile >> d;
		}

		infile.close();
	}
	
	template<typename T>
	void PrintMatrixToConsole(ComputationalMechanicsLibrary::Matrix<T> &_matrix)
	{
		for (size_t i = 0; i < _matrix.Row(); i++)
		{
			for (size_t j = 0; j < _matrix.Column(); j++)
			{
				std::cout << _matrix[i][j];
				if (j != _matrix.Column() - 1)
				{
					std::cout << ',';
				}
			}

			std::cout << std::endl;
		}
	}
}
