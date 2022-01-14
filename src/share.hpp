#ifndef SHARE_H
#define SHARE_H

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/blas.hpp>
#include <bits/stdc++.h>

namespace ublas = boost::numeric::ublas;
typedef int32_t INT;
typedef ublas::row_major ORI;

/**
 * S = M * rho
 * Now each row of S belongs to one party.
 * Each party will have multiple rows.
 * So we create the map party# : list of shares
 * 
 * Rows of S can be grouped into nCt groups of t rows each.
 * The group number is important for reconstruction of s.
 * So each share is represented as {group#, row}
 */
std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> gen_shares(ublas::vector<INT> s, int t, int n, int lk);

/**
 * Given t or more shares, find out the original secret vector
 * Explanation of the `shares` datatype same as in gen_shares
 */
ublas::vector<INT> reconstruct_secret(std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares, int t);

/**
 * Given linear combination of components of t or more shares,
 * find out the result of the same linear combination
 * applied to the original secret vector
 * Explanation of the `pdts` datatype same as in gen_shares.
 * Only difference is that instead of vector<INT>
 * we just have INT which is the linear combination.
 */
INT recostruct_combination(std::map<int, std::vector<std::pair<int, INT>>> pdts, int t);

/**
 * For all the vectors in all the shares,
 * calculate the inner products with a
 */
std::map<int, std::vector<std::pair<int, INT>>> apply_product(std::map<int, std::vector<std::pair<int, ublas::vector<INT>>>> shares, ublas::vector<INT> a);

#endif