//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#include "ribbon_demo.h"

using namespace RIBBON_COMPILE_BLOOM;

void testISoln() {
    using TypeParam = Settings_Coeff128_Homog;

    IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
    IMPORT_RIBBON_IMPL_TYPES(TypeParam);

    Index num_slots = 131072;
    uint32_t num_to_add = 109306;
    double overhead_ratio = 1.0 * num_slots / num_to_add;


    uint32_t max_ibytes = static_cast<uint32_t>(sizeof(ResultRow) * num_slots);

    std::unique_ptr<char[]> idata(new char[max_ibytes]);
    InterleavedSoln isoln(idata.get(), max_ibytes);

    Banding banding;
    Slice slice1("liyou");
    Slice slice2("zhangxin");
    Slice slice3("jiawei");
    Slice slice4("liuliu");
    Slice slice5("liuke");

    banding.SetOrdinalSeed(0);
    banding.Reset(num_slots);
    banding.Add(slice1);
    banding.Add(slice2);
    banding.Add(slice5);
//    banding.Add(slice3);

//    banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end);
//    banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8));

    isoln.BackSubstFrom(banding);

    Hasher hasher;
    bool b1 = isoln.FilterQuery(slice1, hasher);
    bool b2 = isoln.FilterQuery(slice2, hasher);
    bool b3 = isoln.FilterQuery(slice3, hasher);
    bool b4 = isoln.FilterQuery(slice4, hasher);
    bool b5 = isoln.FilterQuery(slice5, hasher);

    bool b00 = false;
}

int main() {
    testISoln();
}
