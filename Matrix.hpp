
	/*	Boost_Competency_Test --- Aayush Surana GSoC '18  */



#include <bits/stdc++.h>

using namespace std;

template <typename T>
class MatrixExp {

	public:

		auto operator[](size_t i) const { return static_cast<T const&>(*this)[i]; }

		size_t size()               const { return static_cast<T const&>(*this).size(); }

		operator T& () 			  { return static_cast<T&>(*this); }

		operator const T& () 	    const { return static_cast<const T&>(*this); }

};

template <typename R>
class Matrix : public MatrixExp< Matrix<R> > { 
	
	valarray<R> data;

	public:

	R operator[](size_t i) const { return data[i]; }
	R &operator[](size_t i)      { return data[i]; }
	size_t size()          const { return data.size(); }

	//data[] contains the matrix elements in the form of 1-D vector
	//data[0] stores the number of rows

	R operator()(size_t row, size_t col) const{
		return data[data[0]*row + col + 1];
	}

	R& operator()(size_t row, size_t col){
		return data[data[0]*row + col + 1];
	}

	Matrix () : data(1)	{data[0]=0;}
//	Matrix (const Matrix<R> &M);
//	Matrix (Matrix<R> &&M);   
	Matrix (size_t n, size_t m) : data(n*m + 1) {assert(n>0 and m>0);	data[0] = n;}

	unsigned int Get_Row_Size();
	unsigned int Get_Column_Size();
	//the addition operation occurs only at the destination when in is reached hence, no memory overhead cost

	template <typename T>
		inline Matrix<R> &operator+=(const MatrixExp<T> &X) {

			//auto addToCurrent = [&]() {
				this->data[0] = X[0];
				for(size_t i = 1; i < X.size(); i++) {
					this->data[i] = this->data[i] + X[i];
				}
				return *this;	
		//	};
		//	return addToCurrent();
		}

	template <typename T>
		inline Matrix<R> &operator=(const MatrixExp<T> &X) {
			
			//auto equalTo = [&]() {
				this->data[0] = X[0];
				for(size_t i = 1; i < X.size(); i++) {
					this->data[i] = X[i];
				}
				return *this;
//			};
//			return equalTo();
		}

	/*template<typename T>
		inline Matrix<R> &operator=(const Matrix<T> &X) {
			this->data[0] = X[0];
			for(size_t i = 1; i < X.size(); i++) {
				this->data[i] = X[i];
			}
			return *this;
		}*/

	template <typename T>
		inline Matrix<R> &operator*=(const MatrixExp<T> &X) {
	  	
			//(total size - 1) divided by the number of rows gives the number of columns of the matrix
		
			assert(((this->data.size()-1)/(this->data[0])) == X[0]);

			valarray<R> ans((X[0])*((X.size()-1)/X[0]) + 1);
			ans[0]=X[0];

				for (size_t i = 0; i < this->data[0]; i++) {
					for (size_t j = 0; j < ((X.size()-1)/(X[0])); j++) {
						for (size_t k = 0; k < (this->data.size()-1)/(this->data[0]); k++) {
	 						ans[X[0]*i + j + 1] += this->data[X[0]*i + k + 1] * X[X[0]*k + j + 1];
	 					}
					}
	 			}

			for (size_t i = 0; i < this->data[0]; i++) {
				for (size_t j = 0; j < ((X.size()-1)/(X[0])); j++) {
					this->data[X[0]*i + j +1] = ans[X[0]*i + j + 1];
				}
			}
			return *this;
		}

	~Matrix()
	{
		data.resize(0);
	}

};
/*
template <typename R>
Matrix<R>::Matrix (const Matrix<R> &M) {
	data.resize( M.data[0]*( ((M.data).size()-1)/M.data[0] ) + 1 );
	for (size_t i = 0; i < M.data[0]; i++) {
		for (size_t j = 0; j < ((M.data).size()-1)/M.data[0]; j++) {
			data[(size_t)(M.data[0])*i + j +1] = M(i,j);
		}
	}
}

template <typename R> 
Matrix<R>::Matrix (Matrix<R> &&M)  {
	data = std::move(M.data);	 
}
*/
template <typename R>
unsigned int Matrix<R>::Get_Row_Size() {
	return data[0];
}

template <typename R>
unsigned int Matrix<R>::Get_Column_Size() {
	return (data.size()-1)/data[0];
}

template <typename E1, typename E2>
class MatrixSum : public MatrixExp< MatrixSum<E1, E2> > {

	E1 const& _u;
	E2 const& _v; 

	public:
	MatrixSum(E1 const& u, E2 const& v) : _u(u), _v(v) {		
		assert(u[0] == v[0] and (u.size()-1) == (v.size() - 1));
	}
	
	auto operator[](size_t i) const { 
		if(i!=0) {
			return _u[i] + _v[i];
		}
		else return _v[i];
	}
	size_t size()	       const { 
		return _v.size(); 
	}
};

template <typename E1, typename E2>
MatrixSum<E1,E2> const
operator+(E1 const& u, E2 const& v) {
	return MatrixSum<E1, E2>(u, v);
}

template <typename E>
Matrix<E> const
operator*(Matrix<E> const& u, Matrix<E> const& v) {
	
	assert(u[0] == (v.size()-1)/v[0]);
	Matrix<E> ans(u[0],((v.size()-1)/v[0]));
	ans[0] = u[0];

	auto multiply = [&]() {
		for (size_t i = 0; i < u[0]; i++) {
			for (size_t j = 0; j < (v.size()-1)/v[0]; j++) {
				for (size_t k = 0; k < (u.size()-1)/u[0]; k++) {
					ans(i,j) += u(i,k)*v(k,j);
				}
			}
		}
		return ans;
	};
	return multiply();
}