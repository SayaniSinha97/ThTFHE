#include "share.hpp"

std::map<std::pair<int, int>, int> ncr_cache;

int ncr(int n, int r)
{
    if (ncr_cache.find({n, r}) == ncr_cache.end()){
        if (r > n || n < 0 || r < 0)
            return 0;
        else{
            if (r == 0 || r == n){
                ncr_cache[{n, r}] = 1;
            }else if (r == 1 || r == n - 1){
                ncr_cache[{n, r}] = n;
            }else{
                ncr_cache[{n, r}] = ncr(n - 1, r) + ncr(n - 1, r - 1);
            }
        }
    }

    return ncr_cache[{n, r}];
}

ublas::matrix<INT, ORI> and_combine(ublas::matrix<INT, ORI>& A, ublas::matrix<INT, ORI>& B)
{
    int rA = A.size1(), cA = A.size2(),
        rB = B.size1(), cB = B.size2();

    int r = rA + rB, c = cA + cB;

    ublas::matrix<INT, ORI> C(r, c);

    C = ublas::zero_matrix<int>(r, c);

    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            if (i < rA){
                if (j == 0 || j == 1){
                    C(i, j) = A(i, 0);
                }else if (j >= 2 && j <= cA){
                    C(i, j) = A(i, j - 1);
                }
            }else{
                if (j == 1){
                    C(i, j) = B(i - rA, 0);
                }else if (j > cA){
                    C(i, j) = B(i - rA, j - cA);
                }
            }
        }
    }

    return C;

}

ublas::matrix<INT, ORI> or_combine(ublas::matrix<INT, ORI>& A, ublas::matrix<INT, ORI>& B)
{
    int rA = A.size1(), cA = A.size2(),
        rB = B.size1(), cB = B.size2();

    int r = rA + rB, c = cA + cB - 1;

    ublas::matrix<INT, ORI> C(r, c);

    C = ublas::zero_matrix<int>(r, c);

    for (int i = 0; i < r; i++){
        for (int j = 0; j < c; j++){
            if (i < rA){
                if (j < cA){
                    C(i, j) = A(i, j);
                }
            }else{
                if (j == 0){
                    C(i, j) = B(i - rA, j);
                }else if (j >= cA){
                    C(i, j) = B(i - rA, j - cA + 1);
                }
            }
        }
    }

    return C;

}

inline void printMatrix(ublas::matrix<INT, ORI> A){
    for (int i = 0; i < A.size1(); i++){
        for (int j = 0; j < A.size2(); j++){
            std::cout << A(i, j) << "\t";
        }
        std::cout << "\n";
    }
}

/**
 * Generate the matrix M according to Benaloh Leichter construction.
 * Here the boolean function is the t-out-of-n majority function.
 * Which is: AND the terms of any t-sized subset of the n variables,
 * then OR all such nCt terms together.
 * For example, the 2-out-of-3 majority function is: x1x2 + x2x3 + x3x1
 */
ublas::matrix<INT, ORI> gen_M(int t, int n){

    int terms = ncr(n, t);
    ublas::matrix<INT, ORI> A(1, 1);
    A(0, 0) = 1;

    std::vector<ublas::matrix<INT, ORI>> res;
    res.push_back(A);
    for (int i = 0; i < t - 1; i++){
        // t terms ANDed
        auto C = and_combine(*(res.rbegin()), A);
        res.push_back(C);
    }
    auto X = *(res.rbegin());
    res.clear();
    res.push_back(X);
    for (int i = 0; i < terms - 1; i++){
        // nCt terms ORed
        auto C = or_combine(*(res.rbegin()), X);
        res.push_back(C);
    }

    auto ans = *(res.rbegin());
    return ans;
}

/**
 * Generates rho for Benaloh Leichter construction.
 * Instead of having one int, our secret is n sized vector
 * So the size of rho is (e, n)
 * The random numbers generated are in [0, 2^lk)
 */
ublas::matrix<INT, ORI> gen_rho(ublas::vector<INT> s, int e, int lk)
{
    int n = s.size();
    ublas::matrix<INT, ORI> rho(e, n);
    rho = ublas::zero_matrix<int>(e, n);
    
    std::default_random_engine gen;
    std::uniform_int_distribution<int> dist(0, (1 << lk));

    for (int i = 0; i < e; i++){
        for (int j = 0; j < n; j++){
            if (i == 0){
                // Copy the secret
                rho(i, j) = s(j);   
            }else{
                // Fill with random
                rho(i, j) = dist(gen);
            }
        }
    }

    return rho;
}

inline int countSetbits(int bs, int n){
    int ans = 0;
    for (int i = 0; i < n; i++){
        ans += bs & 1;
        bs = (bs >> 1);
    }
    return ans;
}


std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> gen_shares(ublas::vector<INT> s, int t, int n, int lk)
{
    auto M = gen_M(t, n);
    int d = M.size1(), e = M.size2();
    int N = s.size();
    auto rho = gen_rho(s, e, lk);

#ifdef SHARE_DEBUG
    printMatrix(M);
    std::cout << "---------\n";
    printMatrix(rho);
    std::cout << "---------\n";
#endif

    static ublas::matrix<INT, ORI> S(d, N);
    S = ublas::zero_matrix<int>(d, N);
    ublas::blas_3::gmm(S, 1, 1, M, rho);

#ifdef SHARE_DEBUG
    printMatrix(S);
    std::cout << "---------\n";
#endif

    static std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares;

    int curr_row = 0;
    for (int bs = 1; bs < (1 << n); bs++){
        if (countSetbits(bs, n) != t) continue;

        int temp = bs;
        for (int i = 0; i < n; i++){
            if ((temp & 1) == 1){
                ublas::vector<INT> row = ublas::row(S, curr_row);
                shares[i].push_back({curr_row, row});
                curr_row++;
            }
            temp = (temp >> 1);
        }

    }

    return shares;

}


ublas::vector<INT> reconstruct_secret(std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares, int t)
{
    int wit = 0;        // How many parties have I witnessed in the map
    std::set<int> groups;
    std::set<int> buffer, buffer2;
    for (auto kv : shares){
        if (wit == t) break;
        
        if (wit == 0){
            for (auto s : kv.second){
                groups.insert(s.first / t);
            }
        }else{
            buffer.clear();
            buffer2.clear();
            for (auto s : kv.second){
                buffer.insert(s.first / t);
            }
            std::set_intersection(groups.begin(), groups.end(),
                                  buffer.begin(), buffer.end(),
                                  std::inserter(buffer2, buffer2.begin()));

            groups.clear();
            for (auto x : buffer2){
                groups.insert(x);
            }
        }

        wit++;
    }

    assert(groups.size() == 1);

    int grp = *(groups.begin());
    wit = 0;
    std::vector<std::pair<int, ublas::vector<INT>>> shareSelect;

    for (auto kv : shares){
        if (wit == t) break;

        for (auto s : kv.second){
            if (s.first / t == grp)
                shareSelect.push_back(s);        
        }

        wit++;
    }

    std::sort(shareSelect.begin(), shareSelect.end(),
        [=](std::pair<int, ublas::vector<INT>>& a, std::pair<int, ublas::vector<INT>>& b){
            return a.first < b.first;
        }
    );

    int N = shareSelect[0].second.size();
    int a = shareSelect.size();

    ublas::matrix<INT, ORI> S(a, N);

    for (int i = 0; i < a; i++){
        for (int j = 0; j < N; j++){
            S(i, j) = shareSelect[i].second(j);
        }
    }

    ublas::matrix<INT, ORI> lambda(1, a);
    lambda(0, 0) = 1;
    for (int i = 1; i < a; i++){
        lambda(0, i) = -1;
    }


    ublas::matrix<INT, ORI> result(1, N);
    result = ublas::zero_matrix<int>(1, N);
    ublas::blas_3::gmm(result, 1, 1, lambda, S);

    static ublas::vector<INT> ans(N);
    for (int i = 0; i < N; i++){
        ans(i) = result(0, i);
    }

#ifdef SHARE_DEBUG
    printMatrix(S);
    std::cout << "--------\n";
    printMatrix(lambda);
    std::cout << "--------\n";
    printMatrix(result);
    std::cout << "--------\n";
#endif

    return ans;
}

/**
 * TODO: This needs a better refactoring. To much copypasta from above
 */
INT recostruct_combination(std::map<int, std::vector<std::pair<int, INT>>> pdts, int t)
{
    int wit = 0;        // How many parties have I witnessed in the map
    std::set<int> groups;
    std::set<int> buffer, buffer2;
    for (auto kv : pdts){
        if (wit == t) break;
        
        if (wit == 0){
            for (auto s : kv.second){
                groups.insert(s.first / t);
            }
        }else{
            buffer.clear();
            buffer2.clear();
            for (auto s : kv.second){
                buffer.insert(s.first / t);
            }
            std::set_intersection(groups.begin(), groups.end(),
                                  buffer.begin(), buffer.end(),
                                  std::inserter(buffer2, buffer2.begin()));

            groups.clear();
            for (auto x : buffer2){
                groups.insert(x);
            }
        }

        wit++;
    }

    assert(groups.size() == 1);

    int grp = *(groups.begin());
    wit = 0;
    std::vector<std::pair<int, INT>> shareSelect;

    for (auto kv : pdts){
        if (wit == t) break;

        for (auto s : kv.second){
            if (s.first / t == grp)
                shareSelect.push_back(s);        
        }

        wit++;
    }

    std::sort(shareSelect.begin(), shareSelect.end(),
        [=](std::pair<int, INT>& a, std::pair<int, INT>& b){
            return a.first < b.first;
        }
    );

    int a = shareSelect.size();

    ublas::matrix<INT, ORI> S(a, 1);

    for (int i = 0; i < a; i++){
        S(i, 0) = shareSelect[i].second;
    }

    ublas::matrix<INT, ORI> lambda(1, a);
    lambda(0, 0) = 1;
    for (int i = 1; i < a; i++){
        lambda(0, i) = -1;
    }


    ublas::matrix<INT, ORI> result(1, 1);
    result = ublas::zero_matrix<int>(1, 1);
    ublas::blas_3::gmm(result, 1, 1, lambda, S);

#ifdef SHARE_DEBUG
    printMatrix(S);
    std::cout << "--------\n";
    printMatrix(lambda);
    std::cout << "--------\n";
    printMatrix(result);
    std::cout << "--------\n";
#endif

    return result(0, 0);
}


std::map<int, std::vector<std::pair<int, INT>>> apply_product(std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares, ublas::vector<INT> a)
{
    std::map<int, std::vector<std::pair<int, INT>>> ans;

    for (auto kv : shares){
        for (auto s : kv.second){
            INT n = ublas::inner_prod(a, s.second);
            ans[kv.first].push_back({s.first, n});
        }
    }

    return ans;
}


inline void printShares(std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares)
{
    for (auto kv : shares){
        std::cout << "Shares for party " << kv.first << "\n";
        for (auto s : kv.second){
            std::cout << "\t" << s.first << " [ ";
            for (int val : s.second){
                std::cout << val << " ";
            }
            std::cout << "]\n";
        }
    }

}

#ifdef SHARE_DEBUG_WITH_MAIN
int main()
{
    int t = 3, n = 5, lk = 3;
    std::vector<int> _s = {1, 0, 1, 1, 0, 0, 1, 1, 0, 1};
    std::vector<int> _a = {5, 3, 2, 1, 4, 5, 6, 7, 9, 10};

    ublas::vector<INT> s(10);
    for (int i = 0; i < s.size(); i++){
        s(i) = _s[i];
    }
    ublas::vector<INT> a(10);
    for (int i = 0; i < a.size(); i++){
        a(i) = _a[i];
    }


    auto shares = gen_shares(s, t, n, lk);
    auto rec_s = reconstruct_secret(shares, t);

    assert(s.size() == rec_s.size());
    
    for (int i = 0; i < s.size(); i++){
        std::cout << s(i) << " ";
    }
    std::cout << "\n";
    for (int i = 0; i < s.size(); i++){
        std::cout << rec_s(i) << " ";
    }
    std::cout << "\n";
    
    for (int i = 0; i < s.size(); i++){
        assert(s(i) == rec_s(i));
    }

    INT actual_pdt = ublas::inner_prod(a, s);
    auto pdts = apply_product(shares, a);
    INT rec_pdt = recostruct_combination(pdts, t);
    std::cout << actual_pdt << " " << rec_pdt << std::endl;
    assert(actual_pdt == rec_pdt);

    return 0;
}
#endif
