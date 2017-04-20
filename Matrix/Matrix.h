#pragma once
#include <cstdlib>
#include <exception>
#include <initializer_list>
#include <vector>

#include <iostream>
#include <fstream>
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
		Matrix(int _row, int _column)
		{
            this->array = new T*[_row];
			this->array[0] = new T[_row*_column]{ 0 };
			for (int i = 1; i < _row; i++)
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
		Matrix(T** arr, int _row, int _column)
		{
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (int j = 0; j < _column; j++)
			{
				this->array[0][j] = arr[0][j];
			}
			for (int i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (int j = 0; j < _column; j++)
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
			int _row = arr.size();
			int _column = arr[0].size();
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (int j = 0; j < _column; j++)
			{
				this->array[0][j] = arr[0][j];
			}
			for (int i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (int j = 0; j < _column; j++)
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
			int _row = arr.size();
			int _column = (arr.begin())->size();
			this->array = new T*[_row];
			this->array[0] = new T[_row*_column];
			for (int j = 0; j < _column; j++)
			{
				this->array[0][j] = *((arr.begin())->begin()+j);
			}
			for (int i = 1; i < _row; i++)
			{
				this->array[i] = this->array[i - 1] + _column;
				for (int j = 0; j < _column; j++)
				{
					this->array[i][j] = *((arr.begin()+i)->begin() + j);
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

        void SetUnit(int _row)
        {
            this->destory();
            this->array = new T*[_row];
            this->array[0] = new T[_row*_row]{0};
            for (int j = 0; j < _row; j++)
            {
                this->array[0][j] = ((T*)arr)[0 * _row + j];
            }
            for(int i=0;i<_row;i++)
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

			for (int i = 0; i < this->row; i++)
			{
				for (int j = 0; j < this->column; j++)
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
		int Row() const
		{
			return this->row;
		}
		/// <summary>
		/// Get the Column.
		/// </summary>
		/// <returns>column</returns>
		int Column() const
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
		T at(int i, int j)
		{
			if (i < 0 || i >= this->row || j < 0 || j >= this->column)
			{
				throw std::exception("out matrix length");
			}

			return this->array[i][j];
		}

		/// <summary>
		/// Determinat
		/// </summary>
		/// <param name="method">The method.</param>
		/// <returns>this matrix's determinat</returns>
		T Determinat(InverseMethod method = InverseMethod::GaussJordan)
		{
			T result;
			try
			{
				switch (method)
				{
				case InverseMethod::GaussJordan:
					result = this->DeterminatGaussJordan();
					break;
				default:
					result = this->DeterminatGaussJordan();
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
			for (int i = 0; i < result.row; i++)
			{
				for (int j = 0; j < result.column; j++)
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
			int n = this->CheckSquare() ? this->row : throw std::exception("not a square");
			//顺便把异常抛了出去
			//制作一个伴随矩阵大小的矩阵
			Matrix<T> result(n, n);

			//存储代数余子式的矩阵（行、列数都比原矩阵少1）
			Matrix<T> temp(n - 1, n - 1);
			//生成伴随矩阵
			for (int i = 0; i < n; i++)
			{
				for (int j = 0; j < n; j++)
				{
					//生成代数余子式
					for (int x = 0; x < n-1; x++)
					{
						for (int y = 0; y < n-1; y++)
						{
							temp.array[x][y] = this->array[x < i ? x : x + 1][y < j ? y : y + 1];
						}
					}
					result.array[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * temp.Determinat();
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
			for (int i = 0; i < result.row; i++)
			{
				for (int j = 0; j < result.column; j++)
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
			for (int i = 0; i < result.row; i++)
			{
				for (int j = 0; j < result.column; j++)
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
			for (int i = 0; i < result.row; i++)
			{
				for (int j = 0; j < result.column; j++)
				{
					for (int k = 0; k < this->column; k++)
					{
						result.array[i][j] += this->array[i][k] * b.array[k][j];
					}
				}
			}
			return result;
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
		Matrix<T>& operator = (const Matrix<T> &b)
		{
			this->FullArray(b);
			return *this;
		}
		/// <summary>
		/// Operator==s the specified b.
		/// </summary>
		/// <param name="b">The b.</param>
		/// <returns></returns>
		bool operator == (const Matrix<T> &b)
		{
			if (this->row != b.row || this->column != b.column)
			{
				if (this->array == nullptr&&b.array == nullptr)
				{
					return true;
				}
				else if (this->array != nullptr&&b.array != nullptr)
				{
					for (int i = 0; i < this->row; i++)
					{
						for (int j = 0; j < this->column; j++)
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
				for (int i = 1; i < n; i++)
				{
					result = result.Multiply(a);
				}

			}
			else if (n == 0)
			{
				result =Matrix<T>(a.Row(), a.Column());
				for (int i = 0; i < result.Row(); i++)
				{
					result.Array()[i][i] = 1;
				}
			}
			else
			{
				Matrix<T> ainv = a.Inverse();
				result = Matrix<T>(ainv);
				for (int i = -1; i > n; i--)
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
			for (int i = 0; i < b.Row(); i++)
			{
				for (int j = 0; j < b.Column(); j++)
				{
					result.Array()[i][j] *= n;
				}
			}
			return result;
		}

		/// <summary>
		/// matrix[i][j]===>>>matrix.array[i][j]
		/// </summary>
		/// <param name="i">The i.</param>
		/// <returns></returns>
		T* operator [] (int i)
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
		int row;
		/// <summary>
		/// The column
		/// </summary>
		int column;
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
			for (int j = 0; j < a.column; j++)
			{
				this->array[0][j] = a.array[0][j];
			}
			for (int i = 1; i < a.row; i++)
			{
				this->array[i] = (this->array[i - 1]) + (a.column);
				for (int j = 0; j < a.column; j++)
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
		void SwapRowSelf(int i, int j)
		{
			T temp;
			for (int k = 0; k < this->column; k++)
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
		void SwapColumnSelf(int i, int j)
		{
			T temp;
			for (int k = 0; k < this->row; k++)
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
			int n = this->CheckSquare() ? this->row : throw std::exception("not a square");
			//顺便把异常抛了出去

			Matrix<T> result(*this);
			int* jr = new int[this->row];
			int* jc = new int[this->column];

			double max;

			for (int k = 0; k < result.row; k++)
			{
				//全选主元
				max = 0;
				jr[k] = k;
				jc[k] = k;

				for (int i = k; i < result.row; i++)
				{
					for (int j = k; j < result.column; j++)
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

				for (int j = 0; j < result.column; j++)
				{
					if (j != k)
					{
						result.array[k][j] *= result.array[k][k];
					}
				}
				for (int i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						for (int j = 0; j < result.column; j++)
						{
							if (j != k)
							{
								result.array[i][j] -= result.array[i][k] * result.array[k][j];
							}
						}
					}
				}
				for (int i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						result.array[i][k] *= -result.array[k][k];
					}
				}

			}

			//恢复
			for (int k = result.row - 1; k >= 0; k--)
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
		/// Determinat by the Gauss Jordan method
		/// </summary>
		/// <returns>
		/// this matrix's determinat
		/// </returns>
		T DeterminatGaussJordan()
		{
			//阶数
			int n = this->CheckSquare() ? this->row : throw std::exception("not a square");			
			//顺便把异常抛了出去
			T resultD = (T)1;


			Matrix<T> result(*this);
			int* jr = new int[this->row];
			int* jc = new int[this->column];

			double max;

			for (int k = 0; k < result.row; k++)
			{
				//全选主元
				max = 0;
				jr[k] = k;
				jc[k] = k;

				for (int i = k; i < result.row; i++)
				{
					for (int j = k; j < result.column; j++)
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

				for (int j = 0; j < result.column; j++)
				{
					if (j != k)
					{
						result.array[k][j] *= result.array[k][k];
					}
				}
				for (int i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						for (int j = 0; j < result.column; j++)
						{
							if (j != k)
							{
								result.array[i][j] -= result.array[i][k] * result.array[k][j];
							}
						}
					}
				}
				for (int i = 0; i < result.row; i++)
				{
					if (i != k)
					{
						result.array[i][k] *= -result.array[k][k];
					}
				}

			}

			//恢复
			for (int k = result.row - 1; k >= 0; k--)
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


	template<typename T>
	void WriteMatrixToCSV(std::string filename, ComputationalMechanicsLibrary::Matrix<T> &matrix)
	{
		std::ofstream outfile(filename);

		for (int i = 0; i < matrix.Row(); i++)
		{
			for (int j = 0; j < matrix.Column(); j++)
			{
				outfile << matrix[i][j] << ",";
			}
			outfile << "\n";
		}

		outfile.close();
	}

	template<typename T>
	void ReadMatrixFromCSV(std::string filename, ComputationalMechanicsLibrary::Matrix<T> &matrix)
	{
		std::ifstream infile(filename);

		char d;
		for (int i = 0; i < matrix.Row(); i++)
		{
			for (int j = 0; j < matrix.Column(); j++)
			{
				infile >> matrix[i][j] >> d;
			}
			//infile >> d;
		}

		infile.close();
	}


}
