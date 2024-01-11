//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#include "ribbon_demo.h"
#include "util/random.h"

using namespace RIBBON_COMPILE_BLOOM;

void testISoln() {
    using TypeParam = Settings_Coeff128_Homog;

    IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
    IMPORT_RIBBON_IMPL_TYPES(TypeParam);
    using KeyGen = TypeParam::KeyGen;

    if (sizeof(CoeffRow) < 8) {
        printf("error in CoeffRow");
        return;
    }

    bool testResult;

    // 不管上面怎么算的，我这边强制设置 大小,这个大小就是hbase 中默认一个chunk的大小
    Index num_slots = 131072;
    uint32_t num_to_add = 109306;
    double overhead_ratio = 1.0 * num_slots / num_to_add;


    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);
    std::string prefix;
    ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

    // Batch that must be added
    std::string added_str = prefix + "added";
    KeyGen keys_begin(added_str, 0);
    KeyGen keys_end(added_str, num_to_add);

    std::string not_str = prefix + "not";
    KeyGen other_keys_begin(not_str, 0);
    KeyGen other_keys_end(not_str, FLAGS_max_check);

    uint32_t max_ibytes = static_cast<uint32_t>(sizeof(ResultRow) * num_slots);

    std::unique_ptr<char[]> idata(new char[max_ibytes]);
    InterleavedSoln isoln(idata.get(), max_ibytes);

    Hasher hasher;
    {
        Banding banding;
        // Traditional solve for a fixed set.
        testResult = (banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));
        if (testResult != true) {
            printf("error in ResetAndFindSeedToSolve\n");
            return;
        } else {
            printf("success in ResetAndFindSeedToSolve and seed = %d\n", banding.GetOrdinalSeed());
        }

        // Also verify that redundant adds are OK (no effect)
        testResult = (
                banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
        if (testResult != true) {
            printf("error in ResetAndFindSeedToSolve\n");
            return;
        } else {
            printf("success in ResetAndFindSeedToSolve which redundant adds and seed = %d\n", banding.GetOrdinalSeed());
        }
        // Now back-substitution
        isoln.BackSubstFrom(banding);
    }

    // Verify keys added
    KeyGen cur = keys_begin;
    while (cur != keys_end) {
        testResult = (isoln.FilterQuery(*cur, hasher));
        if (testResult != true) {
            printf("error in soln query keys added\n");
            return;
        }
        ++cur;
    }
    printf("success in soln query keys added\n");

    // And also check FP rate for isoln
    {
        Index ifp_count = 0;
        cur = other_keys_begin;
        while (cur != other_keys_end) {
            ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
            ++cur;
        }
        double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
        // For expected FP rate, also include false positives due to
        // collisions in Hash value. (Negligible for 64-bit, can matter for
        // 32-bit.)
        double correction =
                FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

        // NOTE: rare violations expected with kHomogeneous
        if (ifp_count >
            FrequentPoissonUpperBound(expected_fp_count + correction)) {
            printf("error in fp up bound islon \n");
        } else {
            printf("success in fp up bound islon\n");
        }

        // FIXME: why sometimes can we slightly "beat the odds"?
        // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
        // (0.95 factor should not be needed)
        if (ifp_count < FrequentPoissonLowerBound(
                0.95 * expected_fp_count + correction)) {
            printf("error in fp low bound islon\n");
        } else {
            printf("success in fp low bound islon\n");
        };
    }
}

int main() {
    testISoln();
}
