//  Copyright (c) Facebook, Inc. and its affiliates. All Rights Reserved.
//  This source code is licensed under both the GPLv2 (found in the
//  COPYING file in the root directory) and Apache 2.0 License
//  (found in the LICENSE.Apache file in the root directory).

#include "ribbon_test.h"
#include "random.h"

template <typename TypesAndSettings>
class RibbonTypeParamTest : public ::testing::Test {};

class RibbonTest : public ::testing::Test {};

using namespace RIBBON_TEST;

using TestTypesAndSettings = ::testing::Types<
    TypesAndSettings_Coeff128, TypesAndSettings_Coeff128Smash,
    TypesAndSettings_Coeff64, TypesAndSettings_Coeff64Smash,
    TypesAndSettings_Coeff64Smash0, TypesAndSettings_Coeff128_Homog,
    TypesAndSettings_Coeff128Smash_Homog, TypesAndSettings_Coeff64_Homog,
    TypesAndSettings_Coeff64Smash_Homog, TypesAndSettings_Result16,
    TypesAndSettings_Result32, TypesAndSettings_IndexSizeT,
    TypesAndSettings_Hash32, TypesAndSettings_Hash32_Result16,
    TypesAndSettings_KeyString, TypesAndSettings_Seed8,
    TypesAndSettings_NoAlwaysOne, TypesAndSettings_AllowZeroStarts,
    TypesAndSettings_Seed64, TypesAndSettings_Rehasher,
    TypesAndSettings_Rehasher_Result16, TypesAndSettings_Rehasher_Result32,
    TypesAndSettings_Rehasher_Seed64, TypesAndSettings_Rehasher32,
    TypesAndSettings_Rehasher32_Coeff64, TypesAndSettings_SmallKeyGen,
    TypesAndSettings_Hash32_SmallKeyGen, TypesAndSettings_Coeff32,
    TypesAndSettings_Coeff32Smash, TypesAndSettings_Coeff16,
    TypesAndSettings_Coeff16Smash>;

TYPED_TEST_CASE(RibbonTypeParamTest, TestTypesAndSettings);

TYPED_TEST(RibbonTypeParamTest, CompactnessAndBacktrackAndFpRate) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;
  using ConfigHelper =
      ROCKSDB_NAMESPACE::ribbon::BandingConfigHelper<TypeParam>;

  if (sizeof(CoeffRow) < 8) {
    ROCKSDB_GTEST_BYPASS("Not fully supported");
    return;
  }

  const auto log2_thoroughness =
      static_cast<uint32_t>(ROCKSDB_NAMESPACE::FloorLog2(FLAGS_thoroughness));

  // 下面给出了中小规模filter的添加key数量范围，原理待研究
  // We are going to choose num_to_add using an exponential distribution,
  // so that we have good representation of small-to-medium filters.
  // Here we just pick some reasonable, practical upper bound based on
  // kCoeffBits or option. // log() 计算以e为底的对数
  const double log_max_add = std::log(
      FLAGS_max_add > 0 ? FLAGS_max_add
                        : static_cast<uint32_t>(kCoeffBits * kCoeffBits) *
                              std::max(FLAGS_thoroughness, uint32_t{32}));

  // This needs to be enough below the minimum number of slots to get a
  // reasonable number of samples with the minimum number of slots.
  const double log_min_add = std::log(0.66 * SimpleSoln::RoundUpNumSlots(1));

  ASSERT_GT(log_max_add, log_min_add);

  const double diff_log_add = log_max_add - log_min_add;

  for (ConstructionFailureChance cs : TypeParam::FailureChanceToTest()) {
    double expected_reseeds;
    switch (cs) {
      default:
        assert(false);
//        FALLTHROUGH_INTENDED;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn2:
        fprintf(stderr, "== Failure: 50 percent\n");
        expected_reseeds = 1.0;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn20:
        fprintf(stderr, "== Failure: 5 percent\n"); //这里应该是写错了，应该是5%，写的是95%。。
        expected_reseeds = 0.053;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn1000:
        fprintf(stderr, "== Failure: 1/1000\n");
        expected_reseeds = 0.001;
        break;
    }

    uint64_t total_reseeds = 0;
    uint64_t total_singles = 0;
    uint64_t total_single_failures = 0;
    uint64_t total_batch = 0;
    uint64_t total_batch_successes = 0;
    uint64_t total_fp_count = 0;
    uint64_t total_added = 0;
    uint64_t total_expand_trials = 0;
    uint64_t total_expand_failures = 0;
    double total_expand_overhead = 0.0;

    uint64_t soln_query_nanos = 0; // soln查询时间
    uint64_t soln_query_count = 0; // soln查询数量
    uint64_t bloom_query_nanos = 0; // bloom查询的时间
    uint64_t isoln_query_nanos = 0; // isoln查询时间
    uint64_t isoln_query_count = 0; // isoln查询数量

    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);

    for (uint32_t i = 0; i < FLAGS_thoroughness; ++i) { // 不太懂这个迭代次数有啥用
      // We are going to choose num_to_add using an exponential distribution
      // as noted above, but instead of randomly choosing them, we generate
      // samples linearly using the golden ratio, which ensures a nice spread
      // even for a small number of samples, and starting with the minimum
      // number of slots to ensure it is tested. // fmod是浮点数取余
      double log_add =
          std::fmod(0.6180339887498948482 * diff_log_add * i, diff_log_add) +
          log_min_add; // 但是可以看到要添加的key的数量是收到迭代次数影响的
      uint32_t num_to_add = static_cast<uint32_t>(std::exp(log_add)); // e的x方

      // Most of the time, test the Interleaved solution storage, but when
      // we do we have to make num_slots a multiple of kCoeffBits. So
      // sometimes we want to test without that limitation.
      bool test_interleaved = (i % 7) != 6; // 迭代次数的作用1：
                                            // 正常来讲num_slots是kCoeffBits的倍数，但有时候我们希望没有这个限制，
                                            // 在迭代次数除7余6的时候，我们取消这个限制

      // Compute num_slots, and re-adjust num_to_add to get as close as possible
      // to next num_slots, to stress that num_slots in terms of construction
      // success. Ensure at least one iteration:
      Index num_slots = Index{0} - 1; // 得到一个Index类型所能表示的最大值
      --num_to_add;
      for (;;) { // 很迷惑，反正通过一些操作，把上面的num_slots重新给了一个较小的值 具体实现在ribbon_config.cc
        Index next_num_slots = SimpleSoln::RoundUpNumSlots(
            ConfigHelper::GetNumSlots(num_to_add + 1, cs));
        if (test_interleaved) {
          next_num_slots = InterleavedSoln::RoundUpNumSlots(next_num_slots);
          // assert idempotent
          EXPECT_EQ(next_num_slots,
                    InterleavedSoln::RoundUpNumSlots(next_num_slots));
        }
        // assert idempotent with InterleavedSoln::RoundUpNumSlots
        EXPECT_EQ(next_num_slots, SimpleSoln::RoundUpNumSlots(next_num_slots));

        if (next_num_slots > num_slots) {
          break;
        }
        num_slots = next_num_slots;
        ++num_to_add;
      }
      assert(num_slots < Index{0} - 1);

      total_added += num_to_add;

      std::string prefix;
      ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

      // Batch that must be added
      std::string added_str = prefix + "added";
      KeyGen keys_begin(added_str, 0);
      KeyGen keys_end(added_str, num_to_add);

      // A couple more that will probably be added
      KeyGen one_more(prefix + "more", 1);
      KeyGen two_more(prefix + "more", 2);

      // Batch that may or may not be added
      uint32_t batch_size =
          static_cast<uint32_t>(2.0 * std::sqrt(num_slots - num_to_add));
      if (batch_size < 10U) {
        batch_size = 0;
      }
      std::string batch_str = prefix + "batch";
      KeyGen batch_begin(batch_str, 0);
      KeyGen batch_end(batch_str, batch_size);

      // Batch never (successfully) added, but used for querying FP rate
      std::string not_str = prefix + "not";
      KeyGen other_keys_begin(not_str, 0);
      KeyGen other_keys_end(not_str, FLAGS_max_check);

      double overhead_ratio = 1.0 * num_slots / num_to_add;
      if (FLAGS_verbose) {
        fprintf(stderr, "Adding(%s) %u / %u   Overhead: %g   Batch size: %u\n",
                test_interleaved ? "i" : "s", (unsigned)num_to_add,
                (unsigned)num_slots, overhead_ratio, (unsigned)batch_size);
      }

      // Vary bytes for InterleavedSoln to use number of solution columns
      // from 0 to max allowed by ResultRow type (and used by SimpleSoln).
      // Specifically include 0 and max, and otherwise skew toward max.
      uint32_t max_ibytes =
          static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
      size_t ibytes;
      // 迭代轮数的作用2：设置不同大小的InterleavedSoln
      if (i == 0) {
        ibytes = 0;
      } else if (i == 1) {
        ibytes = max_ibytes;
      } else {
        // Skewed 偏向于更大的随机数
        ibytes =
            std::max(rnd.Uniformish(max_ibytes), rnd.Uniformish(max_ibytes));
      }
      std::unique_ptr<char[]> idata(new char[ibytes]);
      InterleavedSoln isoln(idata.get(), ibytes);

      SimpleSoln soln;
      Hasher hasher;
      bool first_single;
      bool second_single;
      bool batch_success;
      {
        Banding banding;
        // Traditional solve for a fixed set.
        ASSERT_TRUE(
            banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));

//        banding.printCoeRow((int )num_slots);
//        banding.printResultRow((int )num_slots);

        Index occupied_count = banding.GetOccupiedCount();
        Index more_added = 0;

        if (TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
            batch_size == 0) {
          // 可以看出来backtracking是用来记录banding失败的
          // Homogeneous not compatible with backtracking because add
          // 在开销小于1.01时，我们不应该要求更多
          // doesn't fail. Small overhead ratio too packed to expect more
          first_single = false;
          second_single = false;
          batch_success = false;
        } else {
          // Now to test backtracking, starting with guaranteed fail. By using
          // the keys that will be used to test FP rate, we are then doing an
          // extra check that after backtracking there are no remnants (e.g. in
          // result side of banding) of these entries.
          KeyGen other_keys_too_big_end = other_keys_begin;
          other_keys_too_big_end += num_to_add; // 一般这个num_to_add是大于nums_slots的一半的，所以这里再添加num_to_add一定会超过filter容量
          banding.EnsureBacktrackSize(std::max(num_to_add, batch_size));
          EXPECT_FALSE(banding.AddRangeOrRollBack(other_keys_begin,
                                                  other_keys_too_big_end));// 由于超过容量所以失败
          EXPECT_EQ(occupied_count, banding.GetOccupiedCount());

          // Check that we still have a good chance of adding a couple more
          // individually
          first_single = banding.Add(*one_more);
          second_single = banding.Add(*two_more);
          more_added += (first_single ? 1 : 0) + (second_single ? 1 : 0);
          total_singles += 2U;
          total_single_failures += 2U - more_added;

          // Or as a batch
          batch_success = banding.AddRangeOrRollBack(batch_begin, batch_end);
          ++total_batch;
          if (batch_success) {
            more_added += batch_size;
            ++total_batch_successes;
          }
          EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);
        }

        // Also verify that redundant adds are OK (no effect)
        ASSERT_TRUE(
            banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
        EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);

        // Now back-substitution
        soln.BackSubstFrom(banding);
        if (test_interleaved) {
          isoln.BackSubstFrom(banding);
        }

        Seed reseeds = banding.GetOrdinalSeed();
        total_reseeds += reseeds;

        EXPECT_LE(reseeds, 8 + log2_thoroughness); // 迷惑，其实干嘛要做log运算我是很迷惑的
        if (reseeds > log2_thoroughness + 1) {
          fprintf(
              stderr, "%s high reseeds at %u, %u/%u: %u\n",
              reseeds > log2_thoroughness + 8 ? "ERROR Extremely" : "Somewhat",
              static_cast<unsigned>(i), static_cast<unsigned>(num_to_add),
              static_cast<unsigned>(num_slots), static_cast<unsigned>(reseeds));
        }

        if (reseeds > 0) {
          // "Expand" test: given a failed construction, how likely is it to
          // pass with same seed and more slots. At each step, we increase
          // enough to ensure there is at least one shift within each coeff
          // block.
          ++total_expand_trials;
          Index expand_count = 0;
          Index ex_slots = num_slots;
          banding.SetOrdinalSeed(0);
          for (;; ++expand_count) {
            ASSERT_LE(expand_count, log2_thoroughness);
            ex_slots += ex_slots / kCoeffBits;
            if (test_interleaved) {
              ex_slots = InterleavedSoln::RoundUpNumSlots(ex_slots);
            }
            banding.Reset(ex_slots);
            bool success = banding.AddRange(keys_begin, keys_end);
            if (success) {
              break;
            }
          }
          total_expand_failures += expand_count;
          total_expand_overhead += 1.0 * (ex_slots - num_slots) / num_slots;
        }

        hasher.SetOrdinalSeed(reseeds);
      }
      // soln and hasher now independent of Banding object

      // Verify keys added
      KeyGen cur = keys_begin;
      while (cur != keys_end) {
        ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
        ++cur;
      }
      // We (maybe) snuck these in!
      if (first_single) {
        ASSERT_TRUE(soln.FilterQuery(*one_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*one_more, hasher));
      }
      if (second_single) {
        ASSERT_TRUE(soln.FilterQuery(*two_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*two_more, hasher));
      }
      if (batch_success) {
        cur = batch_begin;
        while (cur != batch_end) {
          ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
          ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
          ++cur;
        }
      }

      // Check FP rate (depends only on number of result bits == solution
      // columns)
      Index fp_count = 0;
      cur = other_keys_begin;
      {
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          bool fp = soln.FilterQuery(*cur, hasher);
          fp_count += fp ? 1 : 0;
          ++cur;
        }
        soln_query_nanos += timer.ElapsedNanos();
        soln_query_count += FLAGS_max_check;
      }
      {
        double expected_fp_count = soln.ExpectedFpRate() * FLAGS_max_check;
        // For expected FP rate, also include false positives due to collisions
        // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
        double correction =
            FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

        // NOTE: rare violations expected with kHomogeneous
        EXPECT_LE(fp_count,
                  FrequentPoissonUpperBound(expected_fp_count + correction)); // 泊松分布对比误差
        EXPECT_GE(fp_count,
                  FrequentPoissonLowerBound(expected_fp_count + correction));
      }
      total_fp_count += fp_count;

      // And also check FP rate for isoln
      if (test_interleaved) {
        Index ifp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
          ++cur;
        }
        isoln_query_nanos += timer.ElapsedNanos();
        isoln_query_count += FLAGS_max_check;
        {
          double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
          // For expected FP rate, also include false positives due to
          // collisions in Hash value. (Negligible for 64-bit, can matter for
          // 32-bit.)
          double correction =
              FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

          // NOTE: rare violations expected with kHomogeneous
          EXPECT_LE(ifp_count,
                    FrequentPoissonUpperBound(expected_fp_count + correction));

          // FIXME: why sometimes can we slightly "beat the odds"?
          // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
          // (0.95 factor should not be needed)
          EXPECT_GE(ifp_count, FrequentPoissonLowerBound(
                                   0.95 * expected_fp_count + correction));
        }
        // Since the bits used in isoln are a subset of the bits used in soln,
        // it cannot have fewer FPs
        // 从这里看出来 isoln对比soln的优势不在fp，在fp上应该是持平的，所以可能在内存上有优化
        EXPECT_GE(ifp_count, fp_count);
      }

      // And compare to Bloom time, for fun
      if (ibytes >= /* minimum Bloom impl bytes*/ 64) {
        Index bfp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          uint64_t h = hasher.GetHash(*cur);
          uint32_t h1 = ROCKSDB_NAMESPACE::Lower32of64(h);
          uint32_t h2 = sizeof(Hash) >= 8 ? ROCKSDB_NAMESPACE::Upper32of64(h)
                                          : h1 * 0x9e3779b9;
          bfp_count +=
              ROCKSDB_NAMESPACE::FastLocalBloomImpl::HashMayMatch(
                  h1, h2, static_cast<uint32_t>(ibytes), 6, idata.get())
                  ? 1
                  : 0;
          ++cur;
        }
        bloom_query_nanos += timer.ElapsedNanos();
        // ensure bfp_count is used
        ASSERT_LT(bfp_count, FLAGS_max_check);
      }
    }

    // "outside" == key not in original set so either negative or false positive
    fprintf(stderr,
            "Simple      outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * soln_query_nanos / soln_query_count);
    fprintf(stderr,
            "Interleaved outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * isoln_query_nanos / isoln_query_count);
    fprintf(stderr,
            "Bloom       outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * bloom_query_nanos / soln_query_count);

    if (TypeParam::kHomogeneous) {
      EXPECT_EQ(total_reseeds, 0U); // 说明kHomogeneous在banding中不会失败，因此不会改变种子的值
    } else {
      double average_reseeds = 1.0 * total_reseeds / FLAGS_thoroughness;
      fprintf(stderr, "Average re-seeds: %g\n", average_reseeds);
      // Values above were chosen to target around 50% chance of encoding
      // success rate (average of 1.0 re-seeds) or slightly better. But 1.15 is
      // also close enough.
      EXPECT_LE(total_reseeds,
                InfrequentPoissonUpperBound(1.15 * expected_reseeds *
                                            FLAGS_thoroughness));
      // Would use 0.85 here instead of 0.75, but
      // TypesAndSettings_Hash32_SmallKeyGen can "beat the odds" because of
      // sequential keys with a small, cheap hash function. We accept that
      // there are surely inputs that are somewhat bad for this setup, but
      // these somewhat good inputs are probably more likely.
      EXPECT_GE(total_reseeds,
                InfrequentPoissonLowerBound(0.75 * expected_reseeds *
                                            FLAGS_thoroughness));
    }

    if (total_expand_trials > 0) {
      double average_expand_failures =
          1.0 * total_expand_failures / total_expand_trials;
      fprintf(stderr, "Average expand failures, and overhead: %g, %g\n",
              average_expand_failures,
              total_expand_overhead / total_expand_trials);
      // Seems to be a generous allowance
      EXPECT_LE(total_expand_failures,
                InfrequentPoissonUpperBound(1.0 * total_expand_trials));
    } else {
      fprintf(stderr, "Average expand failures: N/A\n");
    }

    if (total_singles > 0) {
      double single_failure_rate = 1.0 * total_single_failures / total_singles;
      fprintf(stderr, "Add'l single, failure rate: %g\n", single_failure_rate);
      // A rough bound (one sided) based on nothing in particular
      double expected_single_failures = 1.0 * total_singles /
                                        (sizeof(CoeffRow) == 16 ? 128
                                         : TypeParam::kUseSmash ? 64
                                                                : 32);
      EXPECT_LE(total_single_failures,
                InfrequentPoissonUpperBound(expected_single_failures));
    }

    if (total_batch > 0) {
      // Counting successes here for Poisson to approximate the Binomial
      // distribution.
      // A rough bound (one sided) based on nothing in particular.
      double expected_batch_successes = 1.0 * total_batch / 2;
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_batch_successes);
      fprintf(stderr, "Add'l batch, success rate: %g (>= %g)\n",
              1.0 * total_batch_successes / total_batch,
              1.0 * lower_bound / total_batch);
      EXPECT_GE(total_batch_successes, lower_bound);
    }

    {
      uint64_t total_checked = uint64_t{FLAGS_max_check} * FLAGS_thoroughness;
      double expected_total_fp_count =
          total_checked * std::pow(0.5, 8U * sizeof(ResultRow));
      // For expected FP rate, also include false positives due to collisions
      // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
      double average_added = 1.0 * total_added / FLAGS_thoroughness;
      expected_total_fp_count +=
          total_checked * ExpectedCollisionFpRate(Hasher(), average_added);

      uint64_t upper_bound =
          InfrequentPoissonUpperBound(expected_total_fp_count);
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_total_fp_count);
      fprintf(stderr, "Average FP rate: %g (~= %g, <= %g, >= %g)\n",
              1.0 * total_fp_count / total_checked,
              expected_total_fp_count / total_checked,
              1.0 * upper_bound / total_checked,
              1.0 * lower_bound / total_checked);
      EXPECT_LE(total_fp_count, upper_bound);
      EXPECT_GE(total_fp_count, lower_bound);
    }
  }
}

TYPED_TEST(RibbonTypeParamTest, Extremes2) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;

  size_t bytes = 128 * 1024;
  std::unique_ptr<char[]> buf(new char[bytes]);
  InterleavedSoln isoln(buf.get(), bytes);
  SimpleSoln soln;
  Hasher hasher;
  Banding banding;

  // ########################################
  // Add zero keys to minimal number of slots
  KeyGen begin_and_end("foo", 123);
  KeyGen begin_and_end2("foo", 125);
  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(
      /*slots*/ kCoeffBits, begin_and_end, begin_and_end2, /*first seed*/ 0,
      /* seed mask*/ 0));

  soln.BackSubstFrom(banding);
  isoln.BackSubstFrom(banding);

  for (int i = 123; i < 130; ++i) {
    KeyGen tmp("foo", i);
    bool b1 = isoln.FilterQuery(*tmp, hasher);
    bool b2 = soln.FilterQuery(*tmp, hasher);
    std::cout << "1234" << std::endl;
  }
}

TYPED_TEST(RibbonTypeParamTest, Extremes) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;

  size_t bytes = 128 * 1024;
  std::unique_ptr<char[]> buf(new char[bytes]);
  InterleavedSoln isoln(buf.get(), bytes);
  SimpleSoln soln;
  Hasher hasher;
  Banding banding;

  // ########################################
  // Add zero keys to minimal number of slots
  KeyGen begin_and_end("foo", 123);
  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(
      /*slots*/ kCoeffBits, begin_and_end, begin_and_end, /*first seed*/ 0,
      /* seed mask*/ 0));

  soln.BackSubstFrom(banding);
  isoln.BackSubstFrom(banding);

  // Because there's plenty of memory, we expect the interleaved solution to
  // use maximum supported columns (same as simple solution)
  ASSERT_EQ(isoln.GetUpperNumColumns(), 8U * sizeof(ResultRow));
  ASSERT_EQ(isoln.GetUpperStartBlock(), 0U);

  // Somewhat oddly, we expect same FP rate as if we had essentially filled
  // up the slots.
  KeyGen other_keys_begin("not", 0);
  KeyGen other_keys_end("not", FLAGS_max_check);

  Index fp_count = 0;
  KeyGen cur = other_keys_begin;
  while (cur != other_keys_end) {
    bool isoln_query_result = isoln.FilterQuery(*cur, hasher);
    bool soln_query_result = soln.FilterQuery(*cur, hasher);
    // Solutions are equivalent
    ASSERT_EQ(isoln_query_result, soln_query_result);
    if (!TypeParam::kHomogeneous) {
      // And in fact we only expect an FP when ResultRow is 0
      // (except Homogeneous)
      ASSERT_EQ(soln_query_result, hasher.GetResultRowFromHash(
                                       hasher.GetHash(*cur)) == ResultRow{0});
    }
    fp_count += soln_query_result ? 1 : 0;
    ++cur;
  }
  {
    ASSERT_EQ(isoln.ExpectedFpRate(), soln.ExpectedFpRate());
    double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
    EXPECT_LE(fp_count, InfrequentPoissonUpperBound(expected_fp_count));
    if (TypeParam::kHomogeneous) {
      // Pseudorandom garbage in Homogeneous filter can "beat the odds" if
      // nothing added
    } else {
      EXPECT_GE(fp_count, InfrequentPoissonLowerBound(expected_fp_count));
    }
  }

  // ######################################################
  // Use zero bytes for interleaved solution (key(s) added)

  // Add one key
  KeyGen key_begin("added", 0);
  KeyGen key_end("added", 1);
  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(
      /*slots*/ kCoeffBits, key_begin, key_end, /*first seed*/ 0,
      /* seed mask*/ 0));

  InterleavedSoln isoln2(nullptr, /*bytes*/ 0);

  isoln2.BackSubstFrom(banding);

  ASSERT_EQ(isoln2.GetUpperNumColumns(), 0U);
  ASSERT_EQ(isoln2.GetUpperStartBlock(), 0U);

  // All queries return true
  ASSERT_TRUE(isoln2.FilterQuery(*other_keys_begin, hasher));
  ASSERT_EQ(isoln2.ExpectedFpRate(), 1.0);
}

TEST(RibbonTest, AllowZeroStarts) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypesAndSettings_AllowZeroStarts);
  IMPORT_RIBBON_IMPL_TYPES(TypesAndSettings_AllowZeroStarts);
  using KeyGen = StandardKeyGen;

  InterleavedSoln isoln(nullptr, /*bytes*/ 0);
  SimpleSoln soln;
  Hasher hasher;
  Banding banding;

  KeyGen begin("foo", 0);
  KeyGen end("foo", 1);
  // Can't add 1 entry
  ASSERT_FALSE(banding.ResetAndFindSeedToSolve(/*slots*/ 0, begin, end));

  KeyGen begin_and_end("foo", 123);
  // Can add 0 entries
  ASSERT_TRUE(banding.ResetAndFindSeedToSolve(/*slots*/ 0, begin_and_end,
                                              begin_and_end));

  Seed reseeds = banding.GetOrdinalSeed();
  ASSERT_EQ(reseeds, 0U);
  hasher.SetOrdinalSeed(reseeds);

  // Can construct 0-slot solutions
  isoln.BackSubstFrom(banding);
  soln.BackSubstFrom(banding);

  // Should always return false
  ASSERT_FALSE(isoln.FilterQuery(*begin, hasher));
  ASSERT_FALSE(soln.FilterQuery(*begin, hasher));

  // And report that in FP rate
  ASSERT_EQ(isoln.ExpectedFpRate(), 0.0);
  ASSERT_EQ(soln.ExpectedFpRate(), 0.0);
}

TEST(RibbonTest, RawAndOrdinalSeeds) {
  StandardHasher<TypesAndSettings_Seed64> hasher64;
  StandardHasher<DefaultTypesAndSettings> hasher64_32;
  StandardHasher<TypesAndSettings_Hash32> hasher32;
  StandardHasher<TypesAndSettings_Seed8> hasher8;

  for (uint32_t limit : {0xffU, 0xffffU}) {
    std::vector<bool> seen(limit + 1);
    for (uint32_t i = 0; i < limit; ++i) {
      hasher64.SetOrdinalSeed(i);
      auto raw64 = hasher64.GetRawSeed();
      hasher32.SetOrdinalSeed(i);
      auto raw32 = hasher32.GetRawSeed();
      hasher8.SetOrdinalSeed(static_cast<uint8_t>(i));
      auto raw8 = hasher8.GetRawSeed();
      {
        hasher64_32.SetOrdinalSeed(i);
        auto raw64_32 = hasher64_32.GetRawSeed();
        ASSERT_EQ(raw64_32, raw32);  // Same size seed
      }
      if (i == 0) {
        // Documented that ordinal seed 0 == raw seed 0
        ASSERT_EQ(raw64, 0U);
        ASSERT_EQ(raw32, 0U);
        ASSERT_EQ(raw8, 0U);
      } else {
        // Extremely likely that upper bits are set
        ASSERT_GT(raw64, raw32);
        ASSERT_GT(raw32, raw8);
      }
      // Hashers agree on lower bits
      ASSERT_EQ(static_cast<uint32_t>(raw64), raw32);
      ASSERT_EQ(static_cast<uint8_t>(raw32), raw8);

      // The translation is one-to-one for this size prefix
      uint32_t v = static_cast<uint32_t>(raw32 & limit);
      ASSERT_EQ(raw64 & limit, v);
      ASSERT_FALSE(seen[v]);
      seen[v] = true;
    }
  }
}

namespace {

struct PhsfInputGen {
  PhsfInputGen(const std::string& prefix, uint64_t id) : id_(id) {
    val_.first = prefix;
    ROCKSDB_NAMESPACE::PutFixed64(&val_.first, /*placeholder*/ 0);
  }

  // Prefix (only one required)
  PhsfInputGen& operator++() {
    ++id_;
    return *this;
  }

  const std::pair<std::string, uint8_t>& operator*() {
    // Use multiplication to mix things up a little in the key
    ROCKSDB_NAMESPACE::EncodeFixed64(&val_.first[val_.first.size() - 8],
                                     id_ * uint64_t{0x1500000001});
    // Occasionally repeat values etc.
    val_.second = static_cast<uint8_t>(id_ * 7 / 8);
    return val_;
  }

  const std::pair<std::string, uint8_t>* operator->() { return &**this; }

  bool operator==(const PhsfInputGen& other) const {
    // Same prefix is assumed
    return id_ == other.id_;
  }
  bool operator!=(const PhsfInputGen& other) const {
    // Same prefix is assumed
    return id_ != other.id_;
  }

  uint64_t id_;
  std::pair<std::string, uint8_t> val_;
};

struct PhsfTypesAndSettings : public DefaultTypesAndSettings {
  static constexpr bool kIsFilter = false;
};
}  // namespace

TEST(RibbonTest, PhsfBasic) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(PhsfTypesAndSettings);
  IMPORT_RIBBON_IMPL_TYPES(PhsfTypesAndSettings);

  Index num_slots = 12800;
  Index num_to_add = static_cast<Index>(num_slots / 1.02);

  PhsfInputGen begin("in", 0);
  PhsfInputGen end("in", num_to_add);

  std::unique_ptr<char[]> idata(new char[/*bytes*/ num_slots]);
  InterleavedSoln isoln(idata.get(), /*bytes*/ num_slots);
  SimpleSoln soln;
  Hasher hasher;

  {
    Banding banding;
    ASSERT_TRUE(banding.ResetAndFindSeedToSolve(num_slots, begin, end));

    soln.BackSubstFrom(banding);
    isoln.BackSubstFrom(banding);

    hasher.SetOrdinalSeed(banding.GetOrdinalSeed());
  }

  for (PhsfInputGen cur = begin; cur != end; ++cur) {
    ASSERT_EQ(cur->second, soln.PhsfQuery(cur->first, hasher));
    ASSERT_EQ(cur->second, isoln.PhsfQuery(cur->first, hasher));
  }
}

// Not a real test, but a tool used to build APIs in ribbon_config.h
TYPED_TEST(RibbonTypeParamTest, FindOccupancy) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;

  if (!FLAGS_find_occ) {
    ROCKSDB_GTEST_BYPASS("Tool disabled during unit test runs");
    return;
  }

  KeyGen cur(std::to_string(testing::UnitTest::GetInstance()->random_seed()),
             0);

  Banding banding;
  Index num_slots = InterleavedSoln::RoundUpNumSlots(FLAGS_find_min_slots);
  Index max_slots = InterleavedSoln::RoundUpNumSlots(FLAGS_find_max_slots);
  while (num_slots <= max_slots) {
    std::map<int32_t, uint32_t> rem_histogram;
    std::map<Index, uint32_t> slot_histogram;
    if (FLAGS_find_slot_occ) {
      for (Index i = 0; i < kCoeffBits; ++i) {
        slot_histogram[i] = 0;
        slot_histogram[num_slots - 1 - i] = 0;
        slot_histogram[num_slots / 2 - kCoeffBits / 2 + i] = 0;
      }
    }
    uint64_t total_added = 0;
    for (uint32_t i = 0; i < FLAGS_find_iters; ++i) {
      banding.Reset(num_slots);
      uint32_t j = 0;
      KeyGen end = cur;
      end += num_slots + num_slots / 10;
      for (; cur != end; ++cur) {
        if (banding.Add(*cur)) {
          ++j;
        } else {
          break;
        }
      }
      total_added += j;
      for (auto& slot : slot_histogram) {
        slot.second += banding.IsOccupied(slot.first);
      }

      int32_t bucket =
          static_cast<int32_t>(num_slots) - static_cast<int32_t>(j);
      rem_histogram[bucket]++;
      if (FLAGS_verbose) {
        fprintf(stderr, "num_slots: %u i: %u / %u avg_overhead: %g\r",
                static_cast<unsigned>(num_slots), static_cast<unsigned>(i),
                static_cast<unsigned>(FLAGS_find_iters),
                1.0 * (i + 1) * num_slots / total_added);
      }
    }
    if (FLAGS_verbose) {
      fprintf(stderr, "\n");
    }

    uint32_t cumulative = 0;

    double p50_rem = 0;
    double p95_rem = 0;
    double p99_9_rem = 0;

    for (auto& h : rem_histogram) {
      double before = 1.0 * cumulative / FLAGS_find_iters;
      double not_after = 1.0 * (cumulative + h.second) / FLAGS_find_iters;
      if (FLAGS_verbose) {
        fprintf(stderr, "overhead: %g before: %g not_after: %g\n",
                1.0 * num_slots / (num_slots - h.first), before, not_after);
      }
      cumulative += h.second;
      if (before < 0.5 && 0.5 <= not_after) {
        // fake it with linear interpolation
        double portion = (0.5 - before) / (not_after - before);
        p50_rem = h.first + portion;
      } else if (before < 0.95 && 0.95 <= not_after) {
        // fake it with linear interpolation
        double portion = (0.95 - before) / (not_after - before);
        p95_rem = h.first + portion;
      } else if (before < 0.999 && 0.999 <= not_after) {
        // fake it with linear interpolation
        double portion = (0.999 - before) / (not_after - before);
        p99_9_rem = h.first + portion;
      }
    }
    for (auto& slot : slot_histogram) {
      fprintf(stderr, "slot[%u] occupied: %g\n", (unsigned)slot.first,
              1.0 * slot.second / FLAGS_find_iters);
    }

    double mean_rem =
        (1.0 * FLAGS_find_iters * num_slots - total_added) / FLAGS_find_iters;
    fprintf(
        stderr,
        "num_slots: %u iters: %u mean_ovr: %g p50_ovr: %g p95_ovr: %g "
        "p99.9_ovr: %g mean_rem: %g p50_rem: %g p95_rem: %g p99.9_rem: %g\n",
        static_cast<unsigned>(num_slots),
        static_cast<unsigned>(FLAGS_find_iters),
        1.0 * num_slots / (num_slots - mean_rem),
        1.0 * num_slots / (num_slots - p50_rem),
        1.0 * num_slots / (num_slots - p95_rem),
        1.0 * num_slots / (num_slots - p99_9_rem), mean_rem, p50_rem, p95_rem,
        p99_9_rem);

    num_slots = std::max(
        num_slots + 1, static_cast<Index>(num_slots * FLAGS_find_next_factor));
    num_slots = InterleavedSoln::RoundUpNumSlots(num_slots);
  }
}

// Not a real test, but a tool to understand Homogeneous Ribbon
// behavior (TODO: configuration APIs & tests)
TYPED_TEST(RibbonTypeParamTest, OptimizeHomogAtScale) {
  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;

  if (!FLAGS_optimize_homog) {
    ROCKSDB_GTEST_BYPASS("Tool disabled during unit test runs");
    return;
  }

  if (!TypeParam::kHomogeneous) {
    ROCKSDB_GTEST_BYPASS("Only for Homogeneous Ribbon");
    return;
  }

  KeyGen cur(std::to_string(testing::UnitTest::GetInstance()->random_seed()),
             0);

  Banding banding;
  Index num_slots = SimpleSoln::RoundUpNumSlots(FLAGS_optimize_homog_slots);
  banding.Reset(num_slots);

  // This and "band_ovr" is the "allocated overhead", or slots over added.
  // It does not take into account FP rates.
  double target_overhead = 1.20;
  uint32_t num_added = 0;

  do {
    do {
      (void)banding.Add(*cur);
      ++cur;
      ++num_added;
    } while (1.0 * num_slots / num_added > target_overhead);

    SimpleSoln soln;
    soln.BackSubstFrom(banding);

    std::array<uint32_t, 8U * sizeof(ResultRow)> fp_counts_by_cols;
    fp_counts_by_cols.fill(0U);
    for (uint32_t i = 0; i < FLAGS_optimize_homog_check; ++i) {
      ResultRow r = soln.PhsfQuery(*cur, banding);
      ++cur;
      for (size_t j = 0; j < fp_counts_by_cols.size(); ++j) {
        if ((r & 1) == 1) {
          break;
        }
        fp_counts_by_cols[j]++;
        r /= 2;
      }
    }
    fprintf(stderr, "band_ovr: %g ", 1.0 * num_slots / num_added);
    for (unsigned j = 0; j < fp_counts_by_cols.size(); ++j) {
      double inv_fp_rate =
          1.0 * FLAGS_optimize_homog_check / fp_counts_by_cols[j];
      double equiv_cols = std::log(inv_fp_rate) * 1.4426950409;
      // Overhead vs. information-theoretic minimum based on observed
      // FP rate (subject to sampling error, especially for low FP rates)
      double actual_overhead =
          1.0 * (j + 1) * num_slots / (equiv_cols * num_added);
      fprintf(stderr, "ovr_%u: %g ", j + 1, actual_overhead);
    }
    fprintf(stderr, "\n");
    target_overhead -= FLAGS_optimize_homog_granularity;
  } while (target_overhead > 1.0);
}


struct tmpSettings : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = true;
  static constexpr bool kUseSmash = true;
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};

TEST(RibbonTypeParamTest, CompactnessAndBacktrackAndFpRate2) {
  using TypeParam = tmpSettings;

  ROCKSDB_NAMESPACE::StopWatchNano timeStub(ROCKSDB_NAMESPACE::SystemClock::Default().get());

  IMPORT_RIBBON_TYPES_AND_SETTINGS(TypeParam);
  IMPORT_RIBBON_IMPL_TYPES(TypeParam);
  using KeyGen = typename TypeParam::KeyGen;
  using ConfigHelper =
      ROCKSDB_NAMESPACE::ribbon::BandingConfigHelper<TypeParam>;

  if (sizeof(CoeffRow) < 8) {
    ROCKSDB_GTEST_BYPASS("Not fully supported");
    return;
  }

  const auto log2_thoroughness =
      static_cast<uint32_t>(ROCKSDB_NAMESPACE::FloorLog2(FLAGS_thoroughness));

  // 下面给出了中小规模filter的添加key数量范围，原理待研究
  // We are going to choose num_to_add using an exponential distribution,
  // so that we have good representation of small-to-medium filters.
  // Here we just pick some reasonable, practical upper bound based on
  // kCoeffBits or option. // log() 计算以e为底的对数
  const double log_max_add = std::log(
      FLAGS_max_add > 0 ? FLAGS_max_add
                        : static_cast<uint32_t>(kCoeffBits * kCoeffBits) *
                              std::max(FLAGS_thoroughness, uint32_t{32}));

  // This needs to be enough below the minimum number of slots to get a
  // reasonable number of samples with the minimum number of slots.
  const double log_min_add = std::log(0.66 * SimpleSoln::RoundUpNumSlots(1));

  ASSERT_GT(log_max_add, log_min_add);

  const double diff_log_add = log_max_add - log_min_add;

  for (ConstructionFailureChance cs : TypeParam::FailureChanceToTest()) {
    double expected_reseeds;
    switch (cs) {
      default:
        assert(false);
        FALLTHROUGH_INTENDED;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn2:
        fprintf(stderr, "== Failure: 50 percent\n");
        expected_reseeds = 1.0;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn20:
        fprintf(stderr, "== Failure: 5 percent\n"); //这里应该是写错了，应该是5%，写的是95%。。
        expected_reseeds = 0.053;
        break;
      case ROCKSDB_NAMESPACE::ribbon::kOneIn1000:
        fprintf(stderr, "== Failure: 1/1000\n");
        expected_reseeds = 0.001;
        break;
    }

    uint64_t total_reseeds = 0;
    uint64_t total_singles = 0;
    uint64_t total_single_failures = 0;
    uint64_t total_batch = 0;
    uint64_t total_batch_successes = 0;
    uint64_t total_fp_count = 0;
    uint64_t total_added = 0;
    uint64_t total_expand_trials = 0;
    uint64_t total_expand_failures = 0;
    double total_expand_overhead = 0.0;

    uint64_t soln_query_nanos = 0; // soln查询时间
    uint64_t soln_query_count = 0; // soln查询数量
    uint64_t bloom_query_nanos = 0; // bloom查询的时间
    uint64_t isoln_query_nanos = 0; // isoln查询时间
    uint64_t isoln_query_count = 0; // isoln查询数量

    // Take different samples if you change thoroughness
    ROCKSDB_NAMESPACE::Random32 rnd(FLAGS_thoroughness);

    for (uint32_t i = 0; i < FLAGS_thoroughness; ++i) { // 不太懂这个迭代次数有啥用

      std::cout << std::endl;
      std::cout << "第" << i << "次迭代：" << std::endl;
      timeStub.Start();

      // We are going to choose num_to_add using an exponential distribution
      // as noted above, but instead of randomly choosing them, we generate
      // samples linearly using the golden ratio, which ensures a nice spread
      // even for a small number of samples, and starting with the minimum
      // number of slots to ensure it is tested. // fmod是浮点数取余
      double log_add =
          std::fmod(0.6180339887498948482 * diff_log_add * i, diff_log_add) +
          log_min_add; // 但是可以看到要添加的key的数量是收到迭代次数影响的
      uint32_t num_to_add = static_cast<uint32_t>(std::exp(log_add)); // e的x方

      // Most of the time, test the Interleaved solution storage, but when
      // we do we have to make num_slots a multiple of kCoeffBits. So
      // sometimes we want to test without that limitation.
      bool test_interleaved = (i % 7) != 6; // 迭代次数的作用1：
                                             // 正常来讲num_slots是kCoeffBits的倍数，但有时候我们希望没有这个限制，
                                             // 在迭代次数除7余6的时候，我们取消这个限制

      // Compute num_slots, and re-adjust num_to_add to get as close as possible
      // to next num_slots, to stress that num_slots in terms of construction
      // success. Ensure at least one iteration:
      Index num_slots = Index{0} - 1; // 得到一个Index类型所能表示的最大值
      --num_to_add;
      for (;;) { // 很迷惑，反正通过一些操作，把上面的num_slots重新给了一个较小的值 具体实现在ribbon_config.cc
        Index next_num_slots = SimpleSoln::RoundUpNumSlots(
            ConfigHelper::GetNumSlots(num_to_add + 1, cs));
        if (test_interleaved) {
          next_num_slots = InterleavedSoln::RoundUpNumSlots(next_num_slots);
          // assert idempotent
          EXPECT_EQ(next_num_slots,
                    InterleavedSoln::RoundUpNumSlots(next_num_slots));
        }
        // assert idempotent with InterleavedSoln::RoundUpNumSlots
        EXPECT_EQ(next_num_slots, SimpleSoln::RoundUpNumSlots(next_num_slots));

        if (next_num_slots > num_slots) {
          break;
        }
        num_slots = next_num_slots;
        ++num_to_add;
      }
      assert(num_slots < Index{0} - 1);

      total_added += num_to_add;

      std::cout << "计算num时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 num_to_add = " << num_to_add << ", num_slots = " << num_slots << std::endl;

      std::string prefix;
      ROCKSDB_NAMESPACE::PutFixed32(&prefix, rnd.Next());

      // Batch that must be added
      std::string added_str = prefix + "added";
      KeyGen keys_begin(added_str, 0);
      KeyGen keys_end(added_str, num_to_add);

      // A couple more that will probably be added
      KeyGen one_more(prefix + "more", 1);
      KeyGen two_more(prefix + "more", 2);

      // Batch that may or may not be added
      uint32_t batch_size =
          static_cast<uint32_t>(2.0 * std::sqrt(num_slots - num_to_add));
      if (batch_size < 10U) {
        batch_size = 0;
      }
      std::string batch_str = prefix + "batch";
      KeyGen batch_begin(batch_str, 0);
      KeyGen batch_end(batch_str, batch_size);

      // Batch never (successfully) added, but used for querying FP rate
      std::string not_str = prefix + "not";
      KeyGen other_keys_begin(not_str, 0);
      KeyGen other_keys_end(not_str, FLAGS_max_check);

      double overhead_ratio = 1.0 * num_slots / num_to_add;
      if (FLAGS_verbose) {
        fprintf(stderr, "Adding(%s) %u / %u   Overhead: %g   Batch size: %u\n",
                test_interleaved ? "i" : "s", (unsigned)num_to_add,
                (unsigned)num_slots, overhead_ratio, (unsigned)batch_size);
      }

      // Vary bytes for InterleavedSoln to use number of solution columns
      // from 0 to max allowed by ResultRow type (and used by SimpleSoln).
      // Specifically include 0 and max, and otherwise skew toward max.
      uint32_t max_ibytes =
          static_cast<uint32_t>(sizeof(ResultRow) * num_slots);
      size_t ibytes;
      // 迭代轮数的作用2：设置不同大小的InterleavedSoln
      if (i == 0) {
        ibytes = 0;
      } else if (i == 1) {
        ibytes = max_ibytes;
      } else {
        // Skewed 偏向于更大的随机数
        ibytes =
            std::max(rnd.Uniformish(max_ibytes), rnd.Uniformish(max_ibytes));
      }
      std::unique_ptr<char[]> idata(new char[ibytes]);
      InterleavedSoln isoln(idata.get(), ibytes);

      std::cout << "生成HBaseKey 定义soln时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 其中 overhead_ratio = " << overhead_ratio << std::endl;

      SimpleSoln soln;
      Hasher hasher;
      bool first_single;
      bool second_single;
      bool batch_success;
      {
        Banding banding;
        // Traditional solve for a fixed set.
        ASSERT_TRUE(
            banding.ResetAndFindSeedToSolve(num_slots, keys_begin, keys_end));

        //        banding.printCoeRow((int )num_slots);
        //        banding.printResultRow((int )num_slots);

        Index occupied_count = banding.GetOccupiedCount();
        Index more_added = 0;

        std::cout << "banding的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms， 此时添加了num_to_add个key：" << num_to_add << "banding占据数量" << occupied_count<< std::endl;

        if (TypeParam::kHomogeneous || overhead_ratio < 1.01 ||
            batch_size == 0) {
          // 可以看出来backtracking是用来记录banding失败的
          // Homogeneous not compatible with backtracking because add
          // 在开销小于1.01时，我们不应该要求更多
          // doesn't fail. Small overhead ratio too packed to expect more
          first_single = false;
          second_single = false;
          batch_success = false;
        } else {
          // Now to test backtracking, starting with guaranteed fail. By using
          // the keys that will be used to test FP rate, we are then doing an
          // extra check that after backtracking there are no remnants (e.g. in
          // result side of banding) of these entries.
          KeyGen other_keys_too_big_end = other_keys_begin;
          other_keys_too_big_end += num_to_add; // 一般这个num_to_add是大于nums_slots的一半的，所以这里再添加num_to_add一定会超过filter容量
          banding.EnsureBacktrackSize(std::max(num_to_add, batch_size));
          EXPECT_FALSE(banding.AddRangeOrRollBack(other_keys_begin,
                                                  other_keys_too_big_end));// 由于超过容量所以失败
          EXPECT_EQ(occupied_count, banding.GetOccupiedCount());

          // Check that we still have a good chance of adding a couple more
          // individually
          first_single = banding.Add(*one_more);
          second_single = banding.Add(*two_more);
          more_added += (first_single ? 1 : 0) + (second_single ? 1 : 0);
          total_singles += 2U;
          total_single_failures += 2U - more_added;

          // Or as a batch
          batch_success = banding.AddRangeOrRollBack(batch_begin, batch_end);
          ++total_batch;
          if (batch_success) {
            more_added += batch_size;
            ++total_batch_successes;
          }
          EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);
        }

        std::cout << "测试backtracking的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否被跳过："<<(DefaultTypesAndSettings::kHomogeneous || overhead_ratio < 1.01 ||
                                                                                                                       batch_size == 0 )<< std::endl;

        // Also verify that redundant adds are OK (no effect)
        ASSERT_TRUE(
            banding.AddRange(keys_begin, KeyGen(added_str, num_to_add / 8)));
        EXPECT_LE(banding.GetOccupiedCount(), occupied_count + more_added);

        std::cout << "测试添加冗余key的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms" << std::endl;

        // Now back-substitution
        soln.BackSubstFrom(banding);
        if (test_interleaved) {
          isoln.BackSubstFrom(banding);
        }

        std::cout << "soln求解的时间：" << timeStub.ElapsedNanos(true)/1000/1000 << "ms，是否跳过isoln："<<!test_interleaved<< std::endl;

        Seed reseeds = banding.GetOrdinalSeed();
        total_reseeds += reseeds;

        EXPECT_LE(reseeds, 8 + log2_thoroughness); // 迷惑，其实干嘛要做log运算我是很迷惑的
        if (reseeds > log2_thoroughness + 1) {
          fprintf(
              stderr, "%s high reseeds at %u, %u/%u: %u\n",
              reseeds > log2_thoroughness + 8 ? "ERROR Extremely" : "Somewhat",
              static_cast<unsigned>(i), static_cast<unsigned>(num_to_add),
              static_cast<unsigned>(num_slots), static_cast<unsigned>(reseeds));
        }

        if (reseeds > 0) {
          // "Expand" test: given a failed construction, how likely is it to
          // pass with same seed and more slots. At each step, we increase
          // enough to ensure there is at least one shift within each coeff
          // block.
          ++total_expand_trials;
          Index expand_count = 0;
          Index ex_slots = num_slots;
          banding.SetOrdinalSeed(0);
          for (;; ++expand_count) {
            ASSERT_LE(expand_count, log2_thoroughness);
            ex_slots += ex_slots / kCoeffBits;
            if (test_interleaved) {
              ex_slots = InterleavedSoln::RoundUpNumSlots(ex_slots);
            }
            banding.Reset(ex_slots);
            bool success = banding.AddRange(keys_begin, keys_end);
            if (success) {
              break;
            }
          }
          total_expand_failures += expand_count;
          total_expand_overhead += 1.0 * (ex_slots - num_slots) / num_slots;
        }

        hasher.SetOrdinalSeed(reseeds);

        std::cout << "测试通过小幅度增加空间，是否可以增加banding的成功率的时间（expand）：" << timeStub.ElapsedNanos(true)/1000/1000
                  << "ms，是否被跳过："<<!(reseeds > 0)
                  << ", seed = "<< reseeds<< std::endl;
      }
      // soln and hasher now independent of Banding object

      // Verify keys added
      KeyGen cur = keys_begin;
      while (cur != keys_end) {
        ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
        ++cur;
      }
      // We (maybe) snuck these in!
      if (first_single) {
        ASSERT_TRUE(soln.FilterQuery(*one_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*one_more, hasher));
      }
      if (second_single) {
        ASSERT_TRUE(soln.FilterQuery(*two_more, hasher));
        ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*two_more, hasher));
      }
      if (batch_success) {
        cur = batch_begin;
        while (cur != batch_end) {
          ASSERT_TRUE(soln.FilterQuery(*cur, hasher));
          ASSERT_TRUE(!test_interleaved || isoln.FilterQuery(*cur, hasher));
          ++cur;
        }
      }

      std::cout << "查询已添加的key时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms，是否跳过isoln："<<!(test_interleaved) << std::endl;

      // Check FP rate (depends only on number of result bits == solution
      // columns)
      Index fp_count = 0;
      cur = other_keys_begin;
      {
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          bool fp = soln.FilterQuery(*cur, hasher);
          fp_count += fp ? 1 : 0;
          ++cur;
        }
        soln_query_nanos += timer.ElapsedNanos();
        soln_query_count += FLAGS_max_check;
      }
      {
        double expected_fp_count = soln.ExpectedFpRate() * FLAGS_max_check;
        // For expected FP rate, also include false positives due to collisions
        // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
        double correction =
            FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

        // NOTE: rare violations expected with kHomogeneous
        EXPECT_LE(fp_count,
                  FrequentPoissonUpperBound(expected_fp_count + correction)); // 泊松分布对比误差
        EXPECT_GE(fp_count,
                  FrequentPoissonLowerBound(expected_fp_count + correction));
      }
      total_fp_count += fp_count;

      std::cout << "查询soln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms"<< std::endl;

      // And also check FP rate for isoln
      if (test_interleaved) {
        Index ifp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          ifp_count += isoln.FilterQuery(*cur, hasher) ? 1 : 0;
          ++cur;
        }
        isoln_query_nanos += timer.ElapsedNanos();
        isoln_query_count += FLAGS_max_check;
        {
          double expected_fp_count = isoln.ExpectedFpRate() * FLAGS_max_check;
          // For expected FP rate, also include false positives due to
          // collisions in Hash value. (Negligible for 64-bit, can matter for
          // 32-bit.)
          double correction =
              FLAGS_max_check * ExpectedCollisionFpRate(hasher, num_to_add);

          // NOTE: rare violations expected with kHomogeneous
          EXPECT_LE(ifp_count,
                    FrequentPoissonUpperBound(expected_fp_count + correction));

          // FIXME: why sometimes can we slightly "beat the odds"?
          // 这个测试有时会不通过，因此引入了一个因子来降低期望的误差数，但我感觉这是个好事情？
          // (0.95 factor should not be needed)
          EXPECT_GE(ifp_count, FrequentPoissonLowerBound(
                                   0.95 * expected_fp_count + correction));
        }
        // Since the bits used in isoln are a subset of the bits used in soln,
        // it cannot have fewer FPs
        // 从这里看出来 isoln对比soln的优势不在fp，在fp上应该是持平的，所以可能在内存上有优化
        EXPECT_GE(ifp_count, fp_count);
      }

      std::cout << "查询isoln 未添加的key，并计算fp 时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 是否跳过"<< !test_interleaved<< std::endl;

      // And compare to Bloom time, for fun
      if (ibytes >= /* minimum Bloom impl bytes*/ 64) {
        Index bfp_count = 0;
        cur = other_keys_begin;
        ROCKSDB_NAMESPACE::StopWatchNano timer(
            ROCKSDB_NAMESPACE::SystemClock::Default().get(), true);
        while (cur != other_keys_end) {
          uint64_t h = hasher.GetHash(*cur);
          uint32_t h1 = ROCKSDB_NAMESPACE::Lower32of64(h);
          uint32_t h2 = sizeof(Hash) >= 8 ? ROCKSDB_NAMESPACE::Upper32of64(h)
                                          : h1 * 0x9e3779b9;
          bfp_count +=
              ROCKSDB_NAMESPACE::FastLocalBloomImpl::HashMayMatch(
                  h1, h2, static_cast<uint32_t>(ibytes), 6, idata.get())
                  ? 1
                  : 0;
          ++cur;
        }
        bloom_query_nanos += timer.ElapsedNanos();
        // ensure bfp_count is used
        ASSERT_LT(bfp_count, FLAGS_max_check);
      }
      std::cout << "对比bloom的时间：" << timeStub.ElapsedNanos(true)/1000/1000
                << "ms， 此时第" << i <<  "次迭代结束"<< std::endl;
    }

    // "outside" == key not in original set so either negative or false positive
    fprintf(stderr,
            "Simple      outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * soln_query_nanos / soln_query_count);
    fprintf(stderr,
            "Interleaved outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * isoln_query_nanos / isoln_query_count);
    fprintf(stderr,
            "Bloom       outside query, hot, incl hashing, ns/key: %g\n",
            1.0 * bloom_query_nanos / soln_query_count);

    if (TypeParam::kHomogeneous) {
      EXPECT_EQ(total_reseeds, 0U); // 说明kHomogeneous在banding中不会失败，因此不会改变种子的值
    } else {
      double average_reseeds = 1.0 * total_reseeds / FLAGS_thoroughness;
      fprintf(stderr, "Average re-seeds: %g\n", average_reseeds);
      // Values above were chosen to target around 50% chance of encoding
      // success rate (average of 1.0 re-seeds) or slightly better. But 1.15 is
      // also close enough.
      EXPECT_LE(total_reseeds,
                InfrequentPoissonUpperBound(1.15 * expected_reseeds *
                                            FLAGS_thoroughness));
      // Would use 0.85 here instead of 0.75, but
      // TypesAndSettings_Hash32_SmallKeyGen can "beat the odds" because of
      // sequential keys with a small, cheap hash function. We accept that
      // there are surely inputs that are somewhat bad for this setup, but
      // these somewhat good inputs are probably more likely.
      EXPECT_GE(total_reseeds,
                InfrequentPoissonLowerBound(0.75 * expected_reseeds *
                                            FLAGS_thoroughness));
    }

    std::cout << "测试 reseed 是否符合预期的时间（测试逻辑不太懂）：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms"<< std::endl;

    if (total_expand_trials > 0) {
      double average_expand_failures =
          1.0 * total_expand_failures / total_expand_trials;
      fprintf(stderr, "Average expand failures, and overhead: %g, %g\n",
              average_expand_failures,
              total_expand_overhead / total_expand_trials);
      // Seems to be a generous allowance
      EXPECT_LE(total_expand_failures,
                InfrequentPoissonUpperBound(1.0 * total_expand_trials));
    } else {
      fprintf(stderr, "Average expand failures: N/A\n");
    }

    if (total_singles > 0) {
      double single_failure_rate = 1.0 * total_single_failures / total_singles;
      fprintf(stderr, "Add'l single, failure rate: %g\n", single_failure_rate);
      // A rough bound (one sided) based on nothing in particular
      double expected_single_failures = 1.0 * total_singles /
                                        (sizeof(CoeffRow) == 16 ? 128
                                         : TypeParam::kUseSmash ? 64
                                                                : 32);
      EXPECT_LE(total_single_failures,
                InfrequentPoissonUpperBound(expected_single_failures));
    }

    if (total_batch > 0) {
      // Counting successes here for Poisson to approximate the Binomial
      // distribution.
      // A rough bound (one sided) based on nothing in particular.
      double expected_batch_successes = 1.0 * total_batch / 2;
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_batch_successes);
      fprintf(stderr, "Add'l batch, success rate: %g (>= %g)\n",
              1.0 * total_batch_successes / total_batch,
              1.0 * lower_bound / total_batch);
      EXPECT_GE(total_batch_successes, lower_bound);
    }

    {
      uint64_t total_checked = uint64_t{FLAGS_max_check} * FLAGS_thoroughness;
      double expected_total_fp_count =
          total_checked * std::pow(0.5, 8U * sizeof(ResultRow));
      // For expected FP rate, also include false positives due to collisions
      // in Hash value. (Negligible for 64-bit, can matter for 32-bit.)
      double average_added = 1.0 * total_added / FLAGS_thoroughness;
      expected_total_fp_count +=
          total_checked * ExpectedCollisionFpRate(Hasher(), average_added);

      uint64_t upper_bound =
          InfrequentPoissonUpperBound(expected_total_fp_count);
      uint64_t lower_bound =
          InfrequentPoissonLowerBound(expected_total_fp_count);
      fprintf(stderr, "Average FP rate: %g (~= %g, <= %g, >= %g)\n",
              1.0 * total_fp_count / total_checked,
              expected_total_fp_count / total_checked,
              1.0 * upper_bound / total_checked,
              1.0 * lower_bound / total_checked);
      EXPECT_LE(total_fp_count, upper_bound);
      EXPECT_GE(total_fp_count, lower_bound);
    }
    std::cout << "计算expand、single、batch、avg fp rate是否符合预期的时间（应该没啥计算量） ：" << timeStub.ElapsedNanos(true)/1000/1000
              << "ms， 此时enum["<< cs << "]结束" << std::endl;
  }
}

int main(int argc, char** argv) {
  ROCKSDB_NAMESPACE::port::InstallStackTraceHandler();
  ::testing::InitGoogleTest(&argc, argv);
#ifdef GFLAGS
  ParseCommandLineFlags(&argc, &argv, true);
#endif  // GFLAGS
  return RUN_ALL_TESTS();
}
