# LatticeAlgorithms
格子に関するアルゴリズムがいくつか詰めてあります。誰かのお役に立てば光栄です。

関数の説明
void print_mat(const std::vector<std::vector<auto>> A)
行列Aを表示する関数です。
You can print a matrix A with this function.

void print_vec(const std::vector<auto> v)
ベクトルvを表示する関数です。
You can print a vector A with this function.

double dot(const std::vector<auto> x, const std::vector<auto> y)
二つのベクトルx, yの内積を計算する関数です。
This function computes the inner product of two vectors x and y.

void SizeReduce(std::vector<std::vector<int>>& b, std::vector<std::vector<double>>& mu, const int i, const int j)
基底行列bを部分サイズ基底簡約を行いbとそのGram-Schmidt直交化係数行列muを更新します。
This function partial size basis reduces a lattice basis matirix b, and updates b and its Gram-Schmidt coeficient matrix mu.

void LLLReduce(std::vector<std::vector<int>>& b, const float d = 0.99)
基底行列bをLLL基底簡約し、bを更新します。
This function LLL-reduces a basis matrix b and updates b.z

std::vector<int> enumerate_svp(std::vector<std::vector<double>> mu, std::vector<double> B, double R, std::vector<std::vector<int>> b)
bを基底行列として持つような格子に於いて、その格子上の最短ベクトルをRを数え上げ上界として数え上げます。
This function enumerates the shortest vector on a lattice L(b) with upper bound R.
