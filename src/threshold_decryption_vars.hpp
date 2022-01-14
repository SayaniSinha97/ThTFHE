#pragma once

#include <map>
#include <tfhe/tfhe.h>
#include <tfhe/tfhe_io.h>
#include <tfhe/lwe-functions.h>
#include <tfhe/numeric_functions.h>
#include <tfhe/tlwe_functions.h>

std::map<std::pair<int, int>, int> ncr_cacheT;	/* Stores <<n, r>: C(n, r)> */
std::map<std::pair<int, int>, TLweKey*> shared_key_repo;	/* Stores <<party_id, group_id>: key_share> */
