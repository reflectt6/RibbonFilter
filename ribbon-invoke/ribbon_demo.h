//
// Created by Rain Night on 2023/12/28.
//

#ifndef ROCKSDB_RIBBON_COMPILE_BLOOM_H
#define ROCKSDB_RIBBON_COMPILE_BLOOM_H

#include "util/ribbon_test.h"

namespace RIBBON_COMPILE_BLOOM {

using namespace RIBBON_TEST;
//using namespace HBASE_BLOOM_NAMESPACE;

using namespace ROCKSDB_NAMESPACE;
using namespace ROCKSDB_NAMESPACE::ribbon;

// HBase Bloom
uint32_t chunkByteSizeHint = 128 * 1024; //"io.storefile.bloom.block.size"
// Maximum number of times a Bloom filter can be "folded" if oversized
uint32_t foldFactor = 7; //"io.storefile.bloom.max.fold"
double errorRate = 0.01; //"io.storefile.bloom.error.rate");

// param for test
uint32_t testRowLength = 32; // "hbase row length for test");

// ribbon
uint32_t FLAGS_thoroughness = 5;
uint32_t FLAGS_max_add = 0;
uint32_t FLAGS_min_check = 4000;
uint32_t FLAGS_max_check = 100000;


struct RibbonCompileBloomDefaultSettings : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = false;
  static constexpr bool kUseSmash = false;
};

struct Settings_Coeff128_Homog : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = true;
  static constexpr bool kUseSmash = false;
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};


struct Settings_Coeff128Smash_Homog : public Settings_Coeff128_Homog {
  static constexpr bool kUseSmash = true;
};

struct Settings_Coeff128_Homog2 : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = true;
  static constexpr bool kUseSmash = false;
  using Seed = uint32_t;
  using Key = uint32_t;
  using KeyGen = Hash32KeyGenWrapper<StandardKeyGen>;
  static Hash HashFn(const uint32_t& key, uint32_t raw_seed) {
    // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
    // not pass SmallKeyGen tests below without some seed premixing from
    // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
    return ROCKSDB_NAMESPACE::Hash(
        reinterpret_cast<const char*>(&key), 1, raw_seed);
  }
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};


} // namespace

#endif  // ROCKSDB_RIBBON_COMPILE_BLOOM_H
