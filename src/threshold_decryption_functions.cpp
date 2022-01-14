#include "threshold_decryption_functions.hpp"
#include "threshold_decryption_vars.hpp"

int ncrT(int n, int r){
    if (ncr_cacheT.find({n, r}) == ncr_cacheT.end()){
        if (r > n || n < 0 || r < 0)
            return 0;
        else{
            if (r == 0 || r == n){
                ncr_cacheT[{n, r}] = 1;
            }else if (r == 1 || r == n - 1){
                ncr_cacheT[{n, r}] = n;
            }else{
                ncr_cacheT[{n, r}] = ncrT(n - 1, r) + ncrT(n - 1, r - 1);
            }
        }
    }
    return ncr_cacheT[{n, r}];
}

/* B is a distribution matrix for any single variable x1. So B is a identity matrix of dimension k. A is a distribution matrix of form x1x2...xi. This method returns distribution matrix for x1x2...x(i+1)*/
ublas::matrix<int> andCombineT(ublas::matrix<int>& A, ublas::matrix<int>& B, int k){
	int rA = A.size1(); int cA = A.size2();
	int rB = B.size1(); int cB = B.size2();
	int r = rA + rB; int c = cA + cB;
	ublas::matrix<int> C;
	C = ublas::zero_matrix<int>(r, c);

	#pragma omp parallel for collapse(2)
	for(int col = 0; col < k; col++){
		for(int row = 0; row < rA; row++){
			C(row, col) = A(row, col);
			C(row, col + k) = A(row, col);
		}
	}
	#pragma omp barrier

	#pragma omp parallel for collapse(2)
	for(int col = k; col < 2*k; col++){
		for(int row = rA; row < r; row++){
			C(row, col) = B(row - rA, col - k);
		}
	}
	#pragma omp barrier

	#pragma omp parallel for collapse(2)
	for(int col = 2*k; col < c; col++){
		for(int row = 0; row < rA; row++){
			C(row, col) = A(row, col - k);
		}
	}
	#pragma omp barrier

	A.resize(r, c);
	A = C;
	return A;
}

/* B is distribution matrix of t-sized AND form x1x2...xt. A is distribution matrix of OR-ing i number of such x1x2...xt terms. This method returns distribution matrix of OR-ing (i+1) terms of the form x1x2...xt. */
ublas::matrix<int> orCombineT(ublas::matrix<int>& A, ublas::matrix<int>& B, int k){
	int rA = A.size1(); int cA = A.size2();
	int rB = B.size1(); int cB = B.size2();
	int r = rA + rB; int c = cA + cB - k;
	ublas::matrix<int> C;
	C = ublas::zero_matrix<int>(r, c);

	#pragma omp parallel for collapse(2)
	for(int col = 0; col < k; col++){
		for(int row = 0; row < r; row++){
			if (row < rA){
				C(row, col) = A(row, col);
			}else{
				C(row, col) = B(row - rA, col);
			}
		}
	}
	#pragma omp barrier

	#pragma omp parallel for collapse(2)
	for(int col = k; col < cA; col++){
		for(int row = 0; row < rA; row++){
			C(row, col) = A(row, col);
		}
	}
	#pragma omp barrier

	#pragma omp parallel for collapse(2)
	for(int col = cA; col < c; col++){
		for(int row = rA; row < r; row++){
			C(row, col) = B(row - rA, col - cA + k);
		}
	}
	#pragma omp barrier


	A.resize(r, c);
	A = C;
	return A;
}

// Copies the entirety of src matrix inside dst
// starting from the index (dstR, dstC) in dst
void matrixCopy(ublas::matrix<int> &dst, ublas::matrix<int> &src, int dstR, int dstC){
	// for (int r = 0; r < src.size1(); r++){
	// 	for (int c = 0; c < src.size2(); c++){
	// 		dst(dstR + r, dstC + c) = src(r, c);
	// 	}
	// }
	ublas::matrix_range<ublas::matrix<int>> dstSlice (dst, ublas::range(dstR, dstR + src.size1()), ublas::range(dstC, dstC + src.size2()));
	dstSlice = src;
}

ublas::matrix<int> optAndCombineT(int t, int k){
	ublas::matrix<int> I;
	I = ublas::identity_matrix<int>(k);

	ublas::matrix<int> Mf;
	int kt = k*t;
	Mf = ublas::zero_matrix<int>(kt, kt);

	// Divide the matrix into t*t chunks of k*k sub-matrices
	for (int r = 0; r < t; r++){
		for (int c = 0; c < t; c++){
			if (r == 0 || c == t - r){
				matrixCopy(Mf, I, r*k, c*k);
			}
		}
	}

	return Mf;
}

ublas::matrix<int> optOrCombineT(int k, int t, int l, ublas::matrix<int> &A){
	ublas::matrix<int> F, R;
	F = ublas::zero_matrix<int>(A.size1(), k);
	R = ublas::zero_matrix<int>(A.size1(), A.size2() - k);

	for (int r = 0; r < A.size1(); r++){
		for (int c = 0; c < A.size2(); c++){
			if (c < k){
				F(r, c) = A(r, c);
			}else{
				R(r, c - k) = A(r, c);
			}
		}
	}

	ublas::matrix<int> M;
	M = ublas::zero_matrix<int>(l*k*t, k*(t-1) * l + k);
	for (int i = 0; i < l; i++){
		matrixCopy(M, F, i*k*t, 0);
		matrixCopy(M, R, i*k*t, k + i*k*(t-1));
	}

	return M;
}

/* Build Distribution matrix for OR-ing C(p,t) number of x1x2...xt like terms. */
void buildDistributionMatrix(int t, int k, int p, ublas::matrix<int>& M){
	// ublas::matrix<int> M2;
	// M2 = ublas::identity_matrix<int>(k);
	ublas::matrix<int> M1;
	// M1 = M2;
	// for(int i = 2; i <= t; i++){
	// 	M1 = andCombineT(M1, M2, k);
	// }
	M1 = optAndCombineT(t, k);
	M = optOrCombineT(k, t, ncrT(p, t), M1);
	// for(int i = 2; i <= ncrT(p, t); i++){
	// 	M = orCombineT(M, M1, k);
	// }
}

/* rho is a random binary matrix with first k rows coming from the k rows of the secret key */
ublas::matrix<int> buildRho(int k, int p, int e, TLweKey *key, ublas::matrix<int>& rho){
	int N = key->params->N;
	rho = ublas::matrix<int>(e, N);
	for(int row = 0; row < k; row++){
		for(int col = 0; col < N; col++){
			rho(row,col) = key->key[row].coefs[col];
		}
	}
	std::default_random_engine gen;
    std::uniform_int_distribution<int> dist(0, 1);
	for(int row = k; row < e; row++){
		for(int col = 0; col < N; col++){
			rho(row,col) = dist(gen);
		}
	} 
	return rho;
}

/* Naive Matrix Multiplication */
void multiply(ublas::matrix<int>& C, ublas::matrix<int>& A, ublas::matrix<int>& B){
	int m = A.size1();
	int k = A.size2();
	int n = B.size2();

	double *_A = (double *)malloc(m * k * sizeof(double));
	double *_B = (double *)malloc(k * n * sizeof(double));
	double *_C = (double *)malloc(m * n * sizeof(double));

	for (int r = 0; r < m; r++)
		for (int c = 0; c < k; c++)
			_A[r * k + c] = A(r, c);

	for (int r = 0; r < k; r++)
		for (int c = 0; c < n; c++)
			_B[r * n + c] = B(r, c);

	memset(_C, 0, m * n * sizeof(double));

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1.0, _A, k, _B, n, 0.0, _C, n);

	for (int r = 0; r < m; r++)
		for (int c = 0; c < n; c++)
			C(r, c) = _C[r * n + c];

	free(_A);
	free(_B);
	free(_C);
}

/* Given a group_id, find the party_ids present in (group_id)^th combination out of C(p,t) combinations */
void findParties(std::vector<int>& pt, int gid, int t, int p){
	int mem = 0, tmp;
	pt.clear();
	for(int i = 1; i < p; i++){
		tmp = ncrT(p - i, t - mem -1);
		if(gid > tmp){
			gid -= tmp;
		}
		else{
			pt.push_back(i);
			mem += 1;
		}
		if(mem + (p-i) == t){
			for(int j = i + 1; j <= p; j++){
				pt.push_back(j);
			}
			break;
		}
	}
}

/* Get and store the actual shares of each party */
void distributeShares(ublas::matrix<int>& S, int t, int k, int p, TLweParams *params){
	int r = S.size1(), N = params->N;
	int row = 1, group_id, row_count;
	std::vector<int> parties;
	while(row <= r){
		group_id = ceil(row/(floor)(k*t));
		findParties(parties, group_id, t, p);
		for(int it = 1; it <= t; it++){
			row_count = row + (it - 1) * k;
			TLweKey *key_share = new_TLweKey(params);
			for(int i = 0; i < k; i++){
				for(int j = 0; j < N; j++){
					key_share->key[i].coefs[j] = S(row_count + i - 1, j);
				}
			}
			shared_key_repo[{parties[it-1], group_id}] = key_share;
		}
		row += (k*t);
	}
}

/* Preprocess to share the secret key among p parties */
void shareSecret(int t, int p, TLweKey *key, TLweParams *params){
	int k = key->params->k;
	int N = key->params->N;

	ublas::matrix<int> M;
	buildDistributionMatrix(t, k, p, M);
	int d = M.size1();
	int e = M.size2();

	ublas::matrix<int> rho;
	rho = buildRho(k, p, e, key, rho);

	ublas::matrix<int> shares(d, N);
	multiply(shares, M, rho);	/* shares = M . rho */

	distributeShares(shares, t, k, p, params);
}

/* Given a t-sized list of party-ids compute its rank among total C(p,t) combinations */
int findGroupId(std::vector<int> parties, int t, int p){
	int mem = 0;
	int group_count = 1;
	for(int i = 1; i <= p; i++){
		if(std::find(parties.begin(), parties.end(), i) != parties.end()){
			mem += 1;
		}
		else{
			group_count += ncrT(p - i, t - mem - 1);
		}
		if(mem == t){
			break;
		}
	}
	return group_count;
}

void nonFFTmul(TorusPolynomial *ans, IntPolynomial *S, TorusPolynomial *A, int mod)
{
    TorusPolynomial *temp = new_TorusPolynomial(S->N + A->N);
    for (int i = 0; i < S->N + A->N; i++)
        temp->coefsT[i] = 0;

    for (int i = 0; i < S->N; i++){
        for (int j = 0; j < A->N; j++){
            temp->coefsT[i + j] += ((int64_t)S->coefs[i] * (int64_t)A->coefsT[j]); 
        }
    }

    for (int i = mod; i < S->N + A->N; i++){
        temp->coefsT[i % mod] -= temp->coefsT[i];
    }

    for (int i = 0; i < mod; i++){
        ans->coefsT[i] = temp->coefsT[i];
    }
}

void nonFFTmul2(TorusPolynomial *ans, IntPolynomial *S, TorusPolynomial *A, int mod)
{
	long long *temp = new long long[S->N + A->N];
    for (int i = 0; i < S->N + A->N; i++)
        temp[i] = 0;

    for (int i = 0; i < S->N; i++){
        for (int j = 0; j < A->N; j++){
            temp[i + j] += ((int64_t)S->coefs[i] * (int64_t)A->coefsT[j]); 
        }
    }

    for (int i = mod; i < S->N + A->N; i++){
        temp[i % mod] -= temp[i];
    }

    for (int i = 0; i < mod; i++){
        ans->coefsT[i] = (temp[i] % ((long long)549755809793));
    }

}

void thresholdDecrypt(TorusPolynomial* plaintext, TLweSample *ciphertext, TLweParams* params, std::vector<int> parties, int t, int p, double sd){
	int k = params->k;
	int N = params->N;
	int group_id = findGroupId(parties, t, p);
	TorusPolynomial* acc = new_TorusPolynomial(N);
	for(int i = 0; i < N; i++){
		acc->coefsT[i] = 0;
	}
	#pragma omp parallel for num_threads(t)
	for(int i = 0; i < t; i++){
		TorusPolynomial *tmp = new_TorusPolynomial(N);
		for(int j = 0; j < N; j++){
			tmp->coefsT[j] = 0;
		}
		TorusPolynomial *err = new_TorusPolynomial(N);
		for(int j = 0; j < N; j++){
			err->coefsT[j] = gaussian32(0, sd);
		}
		auto part_key = shared_key_repo[{parties[i], group_id}];
		for(int j = 0; j < k; j++){
			// torusPolynomialAddMulR(tmp, &part_key->key[j], &ciphertext->a[j]);
			nonFFTmul2(tmp, &part_key->key[j], &ciphertext->a[j], ciphertext->a[j].N);
		}
		torusPolynomialAddTo(tmp, err);
		#pragma omp critical
		{
			if(i == 0){
				torusPolynomialAddTo(acc, tmp);
			}
			else{
				torusPolynomialSubTo(acc, tmp);
			}
		}
	}
	#pragma omp barrier
	torusPolynomialSubTo(acc, ciphertext->b);
	int _N = acc->N;
	for (int i = 0; i < _N; i++){
		acc->coefsT[i] = -acc->coefsT[i];
	}

	torusPolynomialCopy(plaintext, acc);
}

void partialDecrypt(TLweSample* ciphertext, TLweParams* params, TorusPolynomial* partial_ciphertext, uint64_t* cycle_counts, int party, std::vector<int> parties, int t, int p, double sd){
	int k = params->k;
	int N = params->N;
	int group_id = findGroupId(parties, t, p);
	unsigned int low, high;
	struct timespec start_time = {0, 0};
    struct timespec end_time = {0, 0};
	auto part_key = shared_key_repo[{parties[party], group_id}];
	clock_gettime(CLOCK_MONOTONIC, &start_time);
	// __asm__ __volatile__("xorl %%eax,%%eax\n cpuid \n" ::: "%eax", "%ebx", "%ecx", "%edx");
 //    __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
	// auto clock_start = (static_cast<uint64_t>(high) << 32) | low;
	TorusPolynomial *partial_cipher = new_TorusPolynomial(N);
	for(int j = 0; j < N; j++){
		partial_cipher->coefsT[j] = 0;
	}
	TorusPolynomial *smudging_err = new_TorusPolynomial(N);
	for(int j = 0; j < N; j++){
		smudging_err->coefsT[j] = gaussian32(0, sd);
	}
	for(int j = 0; j < k; j++){
		torusPolynomialAddMulR(partial_cipher, &part_key->key[j], &ciphertext->a[j]);
	}
	torusPolynomialAddTo(partial_cipher, smudging_err);
	// __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
 //    auto clock_stop = (static_cast<uint64_t>(high) << 32) | low;
 //    cycle_counts[party] += (clock_stop - clock_start);
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	for(int j = 0; j < N; j++){
		partial_ciphertext->coefsT[j] = partial_cipher->coefsT[j];
	}
	// std::cout << "Partial decryption done by party " << parties[party] << ". Cycles taken: " << cycle_counts[party] << std::endl;
	std::cout << "Partial decryption done by party " << parties[party] << ". Time taken: " << (((double)end_time.tv_nsec + 1.0e+9 * end_time.tv_sec) - ((double)start_time.tv_nsec + 1.0e+9 * start_time.tv_sec)) * 1.0e-9 << "s\n";
}


void finalDecrypt(TLweSample* ciphertext, TorusPolynomial** partial_ciphertexts, TLweParams* params, uint64_t* cycle_counts, int party, std::vector<int> parties, int t, int p){
	int N = params->N;
	int result_msg = 0;
	unsigned int low, high;
	struct timespec start_time = {0, 0};
    struct timespec end_time = {0, 0};
    clock_gettime(CLOCK_MONOTONIC, &start_time);
	// __asm__ __volatile__("xorl %%eax,%%eax\n cpuid \n" ::: "%eax", "%ebx", "%ecx", "%edx");
 //    __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
	// auto clock_start = (static_cast<uint64_t>(high) << 32) | low;
	TorusPolynomial* result = new_TorusPolynomial(N);
	torusPolynomialCopy(result, ciphertext->b);
	for(int i = 0; i < t; i++){
		if(i == 0){
			torusPolynomialSubTo(result, partial_ciphertexts[i]);
		}
		else{
			torusPolynomialAddTo(result, partial_ciphertexts[i]);
		}
	}
	for (int i = 0; i < 32; i++){
		result_msg += (modSwitchFromTorus32(result->coefsT[i], MSIZE) << i);
	}
	// __asm__ __volatile__ ("rdtsc" : "=a" (low), "=d" (high));
 //    auto clock_stop = (static_cast<uint64_t>(high) << 32) | low;
 //    cycle_counts[party] += (clock_stop - clock_start);
	clock_gettime(CLOCK_MONOTONIC, &end_time);
	// std::cout << t << "-out-of-" << p << " threshold decryption result by party " << parties[party]  << ": " << result_msg << " Cycle count: " << cycle_counts[party] << std::endl;
	std::cout << t << "-out-of-" << p << " threshold decryption result by party " << parties[party]  << ": " << result_msg << " Time: " << (((double)end_time.tv_nsec + 1.0e+9 * end_time.tv_sec) - ((double)start_time.tv_nsec + 1.0e+9 * start_time.tv_sec)) * 1.0e-9 << "s\n";
}