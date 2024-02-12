#include<iostream>
#include<vector>
#include<cmath>
#include<random>
#include<tuple>
#pragma GCC target("avx2")
#define EPSILON 0.0000000000000000000000000000001

//==================
//Prints a matrix A.
//==================
void print_mat(const std::vector<std::vector<auto>> A){
  puts("[");
  for(long i = 0; i < A.size(); ++i){
    printf("[");
    for(long j = 0; j < A.at(0).size(); ++j) std::cout << A.at(i).at(j) << " ";
    puts("]");
  }
  puts("]");
}

//==================
//Prints a vector v.
//==================
void print_vec(const std::vector<auto> v){
  printf("[");
  for(long i = 0; i < v.size(); ++i) std::cout << v.at(i) << " ";
  puts("]");
}

//===============
//Tests if a = b.
//===============
bool equal_RR(auto a, auto b){
	if(std::abs(a - b) < EPSILON) return true;
	return false;
}

//==============================
//Generates zero integer vector.
//==============================
std::vector<int> zero_vector_ZZ(const int n){const std::vector<int> z(n); return z;}

//=====================================
//Tests if the vector v is zero vector.
//=====================================
bool IsZero(const std::vector<auto> v){
	const int n = v.size();
	for(int i = 0; i < n; ++i)
		if(std::abs(v.at(i)) > EPSILON) return false;
	return true;
}

//=========================================
//Computes inner product of vectors x and y
//=========================================
double dot(const std::vector<auto> x, const std::vector<auto> y){
  double z = 0.0;
  const int n = x.size();
  for(int i = 0; i < n; ++i) z += x.at(i) * y.at(i);
  return z;
}

//============================================
//Computes a product of a vector and a matrix.
//============================================
std::vector<int> MulVecMat(const std::vector<auto> x, const std::vector<std::vector<auto>> A){
  const int m = A.at(0).size(), n = x.size(); int i, j;
  std::vector<int> v(m);
  for(i = 0; i < m; ++i) for(j = 0; j < n; ++j) v.at(i) += x.at(j) * A.at(j).at(i);
  return v;
}

std::vector<int> CoefToLattice(const std::vector<auto> v, const std::vector<std::vector<auto>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j;
  std::vector<int> x(m);
  for(i = 0; i < n; ++i){
    for(j = 0; j < m; ++j) x.at(j) += v.at(i) * b.at(i).at(j);
  }
  return x;
}

//==========================================
//Gram_Schmidt's orthogonalization algorithm
//==========================================
std::tuple<std::vector<double>, std::vector<std::vector<double>>> Gram_Schmidt_squared(const std::vector<std::vector<int>> b){
  const int n = b.size(), m = b.at(0).size(); int i, j, k;
  std::vector<double> B(n);
  std::vector<std::vector<double>> GSOb(n, std::vector<double>(m)), mu(n, std::vector<double>(n));
  for(i = 0; i < n; ++i){
  	mu.at(i).at(i) = 1;//identity matrix
    for(j = 0; j < m; ++j) GSOb.at(i).at(j) = b.at(i).at(j);
    for(j = 0; j < i; ++j){
      mu.at(i).at(j) = dot(b.at(i), GSOb.at(j)) / dot(GSOb.at(j), GSOb.at(j));
      for(k = 0; k < m; ++k) GSOb.at(i).at(k) -= mu.at(i).at(j) * GSOb.at(j).at(k);
    }
    B.at(i) = dot(GSOb.at(i), GSOb.at(i));
  }
  return std::forward_as_tuple(B, mu);
}

void SizeReduce(std::vector<std::vector<int>>& b, std::vector<std::vector<double>>& mu, const int i, const int j){
  int q;
  const int m = b.at(0).size();
  if(mu.at(i).at(j) > 0.5 || mu.at(i).at(j) < -0.5){
    q = round(mu.at(i).at(j));
    for(int k = 0; k < m; ++k) b.at(i).at(k) -= q * b.at(j).at(k);
    for(int k = 0; k <= j; ++k) mu.at(i).at(k) -= mu.at(j).at(k) * q;
  }
}

void LLLReduce(std::vector<std::vector<int>>& b, const float d = 0.99){
  const int n = b.size(), m = b.at(0).size(); int j, i;
  std::vector<std::vector<double>> mu;
  std::vector<double> B; std::tie(B, mu) = Gram_Schmidt_squared(b);
  int tmp;
  for(int k = 1; k < n;){
    for(j = k - 1; j > -1; --j) SizeReduce(b, mu, k, j);

    //Checks if the lattice basis matrix b satisfies Lovasz condition.
    if(B.at(k) >= (d - mu.at(k).at(k - 1) * mu.at(k).at(k - 1)) * B.at(k - 1)) ++k;
    else{
      for(i = 0; i < m; ++i){tmp = b.at(k - 1).at(i); b.at(k - 1).at(i) = b.at(k).at(i); b.at(k).at(i) = tmp;}
      std::tie(B, mu) = Gram_Schmidt_squared(b);
      k = std::max(k - 1, 1);
    }
  }
}

//=======================================================
//Enumerate a lattice vector whose norm is less than R^2.
//=======================================================
std::vector<int> ENUM(std::vector<std::vector<double>> mu, std::vector<double> B, double R){
	const int n = B.size(), n1 = n + 1;
	int last_nonzero = 0, k = 0, i;
	double tmp;
	std::vector<int> r(n1), v(n), w(n);
	std::vector<double> c(n), rho(n1);
	std::vector<std::vector<double>> sigma(n1, std::vector<double>(n));

	v.at(0) = 1;//Avoids zero vector.
	for(i = 0; i < n; ++i) r.at(i) = i;
	while(true){
		tmp = v.at(k) - c.at(k);
		rho.at(k) = rho.at(k + 1) + tmp * tmp * B.at(k);
		if(rho.at(k) <= R){
			if(k == 0) return v;
			--k;
			//std::cout << k << std::endl;
			r.at(k) = std::max(r.at(k), r.at(k + 1));
			for(i = r.at(k); i > k; --i) sigma.at(i).at(k) = sigma.at(i + 1).at(k) + mu.at(i).at(k) * v.at(i);
			c.at(k) = -sigma.at(k + 1).at(k);
			v.at(k) = round(c.at(k));
			w.at(k) = 1;
		}else{
			++k;
			if(k == n) return zero_vector_ZZ(n);
			r.at(k) = k;
			if(k >= last_nonzero){
				last_nonzero = k;
				++v.at(k);
			}else{
				if(v.at(k) > c.at(k)) v.at(k) -= w.at(k); else v.at(k) += w.at(k);
				++w.at(k);
			}
		}
	}
}

//===============================
//Enumerates the shortest vector.
//===============================
std::vector<int> enumerate_svp(std::vector<std::vector<double>> mu, std::vector<double> B, double R, std::vector<std::vector<int>> b){
	const int n = B.size();
	std::vector<int> enum_v(n), pre_enum_v;
	while(true){
		pre_enum_v = enum_v;
		enum_v = ENUM(mu, B, R);
		if(IsZero(enum_v)) return CoefToLattice(pre_enum_v, b);
		R = dot(CoefToLattice(enum_v, b), CoefToLattice(enum_v, b)) - 1;
	}
}

int main(int argc, char **argv){
	const int n = atoi(argv[1]);
	int j;
	double R;
	std::vector<double> B;
	std::vector<std::vector<int>> b(n, std::vector<int>(n)), b_copy;
	std::vector<std::vector<double>> mu;
	std::vector<int> svp, v(n);

	//Generates a random lattice basis matrix
	for(int i = 0; i < n; ++i){
		std::random_device rnd;
		std::mt19937 mt(rnd());
		std::uniform_int_distribution<> rand_bound(-99999, 99999); 
		b.at(i).at(i) = 1;
		b.at(i).at(0) = rand_bound(mt);

		v.at(i) = std::abs(rand_bound(mt)) % n;
	}
	//b_copy = b;

	puts("Input lattice basis matrix:");
	print_mat(b);

	//LLLreduce
	LLLReduce(b);

	//Enumerates the shortest vector
	R = dot(b.at(0), b.at(0)) + 1;
	std::tie(B, mu) = Gram_Schmidt_squared(b);
	puts("The shortest vector on the lattice:");
	svp = enumerate_svp(mu, B, R, b);
	print_vec(svp);

	puts("LLL-reduced basis matrix of the input matrix:");
	print_mat(b);
	return 0;
}
