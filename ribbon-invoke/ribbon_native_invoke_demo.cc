//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#include "ribbon_demo.h"
#include "util/random.h"

using namespace RIBBON_COMPILE_BLOOM;

void RepeatCompactnessAndBacktrackAndFpRate2() {
    using TypeParam = Settings_Coeff128_Homog;

//    ROCKSDB_NAMESPACE::StopWatchNano timeStub(
//            ROCKSDB_NAMESPACE::SystemClock::Default().get());

    IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
    IMPORT_RIBBON_IMPL_TYPES(TypeParam);
    using KeyGen = TypeParam::KeyGen;

    if (sizeof(CoeffRow) < 8) {
        printf("error in CoeffRow");
    }

    uint64_t total_reseeds = 0;
    uint64_t total_fp_count = 0;
    uint64_t total_added = 0;

    uint64_t soln_query_nanos = 0;   // soln查询时间
    uint64_t soln_query_count = 0;   // soln查询数量
    uint64_t isoln_query_nanos = 0;  // isoln查询时间
    uint64_t isoln_query_count = 0;  // isoln查询数量

    bool testResult;

    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);

    for (uint32_t i = 0; i < FLAGS_thoroughness; ++i) {  // 不太懂这个迭代次数有啥用

        std::cout << std::endl;
        std::cout << "第" << i << "次迭代：" << std::endl;
//        timeStub.Start();

        // 不管上面怎么算的，我这边强制设置 大小,这个大小就是hbase 中默认一个chunk的大小
        Index num_slots = 131072;
        uint32_t num_to_add = 109306;

        bool test_interleaved = (i % 7) != 6;

        total_added += num_to_add;

        std::string prefix;
        ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

//        std::cout << "计算num时间：" << timeStub.ElapsedNanos(true) / 1000 / 1000
//                  << "ms， 其中 num_to_add = " << num_to_add
//                  << ", num_slots = " << num_slots << std::endl;

        // Batch that must be added
        std::string added_str = prefix + "added";
        KeyGen keys_begin(added_str, 0);
        KeyGen keys_end(added_str, num_to_add);

        std::string not_str = prefix + "not";
        KeyGen other_keys_begin(not_str, 0);
        KeyGen other_keys_end(not_str, FLAGS_max_check);

        double overhead_ratio = 1.0 * num_slots / num_to_add;


        // Vary bytes for InterleavedSoln to use number of solution columns
        // from 0 to max allowed by ResultRow type (and used by SimpleSoln).
        // Specifically include 0 and max, and otherwise skew toward max.
        uint32_t max_ibytes = static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
        size_t ibytes;
        // 迭代轮数的作用2：设置不同大小的InterleavedSoln
        if (i == 0) {
            ibytes = 0;
        } else if (i == 1) {
            ibytes = max_ibytes;
        } else {
            // Skewed 偏向于更大的随机数
            ibytes = std::max(rnd.Uniformish(max_ibytes), rnd.Uniformish(max_ibytes));
        }
        std::unique_ptr<char[]> idata(new char[ibytes]);
        InterleavedSoln isoln(idata.get(), ibytes);

//        std::cout << "生成HBaseKey 定义soln时间："
//                  << timeStub.ElapsedNanos(true) / 1000 / 1000
//                  << "ms， 其中 overhead_ratio = " << overhead_ratio << std::endl;
        SimpleSoln soln;
        Hasher hasher;
        {
            Banding banding;
            // Traditional solve for a fixed set.
            testResult = (banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));
            if (testResult != true) {
                printf("error in ResetAndFindSeedToSolve");
                return;
            }

            //        banding.printCoeRow((int )num_slots);
            //        banding.printResultRow((int )num_slots);

            Index occupied_count = banding.GetOccupiedCount();
            Index more_added = 0;

//            std::cout << "banding的时间："
//                      << timeStub.ElapsedNanos(true) / 1000 / 1000
//                      << "ms， 此时添加了num_to_add个key：" << num_to_add
//                      << "banding占据数量" << occupied_count << std::endl;

            // Also verify that redundant adds are OK (no effect)
            testResult = (
                    banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
            if (testResult != true) {
                printf("error in ResetAndFindSeedToSolve");
                return;
            }
            if (banding.GetOccupiedCount() > (occupied_count + more_added)) {
                printf("error in banding.GetOccupiedCount");
                return;
            }

//            std::cout << "测试添加冗余key的时间："
//                      << timeStub.ElapsedNanos(true) / 1000 / 1000 << "ms"
//                      << std::endl;

            // Now back-substitution
            soln.BackSubstFrom(banding);
            if (test_interleaved) {
                isoln.BackSubstFrom(banding);
            }

//            std::cout << "soln求解的时间："
//                      << timeStub.ElapsedNanos(true) / 1000 / 1000
//                      << "ms，是否跳过isoln：" << !test_interleaved << std::endl;

            Seed reseeds = banding.GetOrdinalSeed();
            total_reseeds += reseeds;


            hasher.SetOrdinalSeed(reseeds);
//            std::cout << "测试通过小幅度增加空间，是否可以增加banding的成功率的时间（"
//                         "expand）："
//                      << timeStub.ElapsedNanos(true) / 1000 / 1000
//                      << "ms，是否被跳过：" << !(reseeds > 0)
//                      << ", seed = " << reseeds << std::endl;
        }
        // soln and hasher now independent of Banding object

        // Verify keys added
        KeyGen cur = keys_begin;
        while (cur != keys_end) {
            testResult = (soln.FilterQuery(*cur, hasher));
            if (testResult != true) {
                printf("error in soln query keys added");
                return;
            }
            testResult = (!test_interleaved || isoln.FilterQuery(*cur, hasher));
            if (testResult != true) {
                printf("error in soln query keys added");
                return;
            }
            ++cur;
        }
//        std::cout << "查询已添加的key时间："
//                  << timeStub.ElapsedNanos(true) / 1000 / 1000
//                  << "ms，是否跳过isoln：" << !(test_interleaved) << std::endl;

        // Check FP rate (depends only on number of result bits == solution
        // columns)
        Index fp_count = 0;
        cur = other_keys_begin;
        {
//            ROCKSDB_NAMESPACE::StopWatchNano timer(
//                    ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
            while (cur != other_keys_end) {
                bool fp = soln.FilterQuery(*cur, hasher);
                fp_count += fp ? 1 : 0;
                ++cur;
            }
//            soln_query_nanos += timer.ElapsedNanos();
            soln_query_count += FLAGS_max_check;
        }
        {
            double expected_fp_count = soln.ExpectedFpRate() * FLAGS_max_check;
            // For expected FP rate, also include false positives due to collisions
            // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
            double correction =
                    FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

            // NOTE: rare violations expected with kHomogeneous
            if (fp_count >
                      FrequentPoissonUpperBound(expected_fp_count +
                                                correction))  // 泊松分布对比误差
            {
                printf("error in fp upper bound");
            }
            if (fp_count <
                      FrequentPoissonLowerBound(expected_fp_count + correction)) {
                printf("error in fp lower bound");
            }
        }
        total_fp_count += fp_count;
//        std::cout << "查询soln 未添加的key，并计算fp 时间："
//                  << timeStub.ElapsedNanos(true) / 1000 / 1000 << "ms" << std::endl;

        // And also check FP rate for isoln
        if (test_interleaved) {
            Index ifp_count = 0;
            cur = other_keys_begin;
//            ROCKSDB_NAMESPACE::StopWatchNano timer(
//                    ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
            while (cur != other_keys_end) {
                ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
                ++cur;
            }
//            isoln_query_nanos += timer.ElapsedNanos();
            isoln_query_count += FLAGS_max_check;
            {
                double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
                // For expected FP rate, also include false positives due to
                // collisions in Hash value. (Negligible for 64-bit, can matter for
                // 32-bit.)
                double correction =
                        FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

                // NOTE: rare violations expected with kHomogeneous
                if (ifp_count >
                          FrequentPoissonUpperBound(expected_fp_count + correction)) {
                    printf("error in fp up bound islon");
                }

                // FIXME: why sometimes can we slightly "beat the odds"?
                // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
                // (0.95 factor should not be needed)
                if (ifp_count < FrequentPoissonLowerBound(
                        0.95 * expected_fp_count + correction)) {
                    printf("error in fp low bound islon");
                }
            }
            // Since the bits used in isoln are a subset of the bits used in soln,
            // it cannot have fewer FPs
            // 从这里看出来
            // isoln对比soln的优势不在fp，在fp上应该是持平的，所以可能在内存上有优化
            if (ifp_count < fp_count) {
                printf("error in fp compare with isoln and soln");
            }
        }
//        std::cout << "查询isoln 未添加的key，并计算fp 时间："
//                  << timeStub.ElapsedNanos(true) / 1000 / 1000 << "ms， 是否跳过"
//                  << !test_interleaved << std::endl;
    }
}

int main() {
    RepeatCompactnessAndBacktrackAndFpRate2();
}
