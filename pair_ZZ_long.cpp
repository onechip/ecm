#include "pair_ZZ_long.h"

NTL_START_IMPL;

NTL_pair_impl(ZZ,long,pair_ZZ_long);
NTL_pair_io_impl(ZZ,long,pair_ZZ_long);
NTL_pair_eq_impl(ZZ,long,pair_ZZ_long);

NTL_vector_impl(pair_ZZ_long,vec_pair_ZZ_long);
NTL_io_vector_impl(pair_ZZ_long,vec_pair_ZZ_long);
NTL_eq_vector_impl(pair_ZZ_long,vec_pair_ZZ_long);

NTL_END_IMPL;
