//
// Created by Rain Night on 2023/12/28.
//

#ifndef ROCKSDB_RIBBON_TEST_H
#define ROCKSDB_RIBBON_TEST_H
#include <stdint.h>
#include "coding.h"
#include "hash.h"
#include "ribbon_config.h"
#include "ribbon_impl.h"

using ROCKSDB_NAMESPACE::ribbon::ExpectedCollisionFpRate;
using ROCKSDB_NAMESPACE::ribbon::StandardHasher;
using ROCKSDB_NAMESPACE::ribbon::StandardRehasherAdapter;

namespace RIBBON_TEST {

// Different ways of generating keys for testing

// Generate semi-sequential keys
struct StandardKeyGen {
  StandardKeyGen(const std::string& prefix, uint64_t id)
      : id_(id), str_(prefix) {
    ROCKSDB_NAMESPACE::PutFixed64(&str_, /*placeholder*/ 0);
  }

  // Prefix (only one required)
  StandardKeyGen& operator++() {
    ++id_;
    return *this;
  }

  StandardKeyGen& operator+=(uint64_t i) {
    id_ += i;
    return *this;
  }

  const std::string& operator*() {
    // Use multiplication to mix things up a little in the key
    ROCKSDB_NAMESPACE::EncodeFixed64(&str_[str_.size() - 8],
                                     id_ * uint64_t{0x1500000001});
    return str_;
  }

  bool operator==(const StandardKeyGen& other) const {
    // Same prefix is assumed
    return id_ == other.id_;
  }
  bool operator!=(const StandardKeyGen& other) const {
    // Same prefix is assumed
    return id_ != other.id_;
  }

  uint64_t id_;
  std::string str_;
};

using ROCKSDB_NAMESPACE::ribbon::ConstructionFailureChance;

const std::vector<ConstructionFailureChance> kFailureOnly50Pct = {
    ROCKSDB_NAMESPACE::ribbon::kOneIn2};

const std::vector<ConstructionFailureChance> kFailureOnlyRare = {
    ROCKSDB_NAMESPACE::ribbon::kOneIn1000};

const std::vector<ConstructionFailureChance> kFailureAll = {
    ROCKSDB_NAMESPACE::ribbon::kOneIn2, ROCKSDB_NAMESPACE::ribbon::kOneIn20,
    ROCKSDB_NAMESPACE::ribbon::kOneIn1000};

struct DefaultTypesAndSettings {
  using CoeffRow = ROCKSDB_NAMESPACE::Unsigned128;
  using ResultRow = uint8_t;
  using Index = uint32_t;
  using Hash = uint64_t;
  using Seed = uint32_t;
  using Key = ROCKSDB_NAMESPACE::Slice;
  static constexpr bool kIsFilter = true;
  static constexpr bool kHomogeneous = false;
  static constexpr bool kFirstCoeffAlwaysOne = true;
  static constexpr bool kUseSmash = false;
  static constexpr bool kAllowZeroStarts = false;
  static Hash HashFn(const Key& key, uint64_t raw_seed) {
    // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
    // not pass SmallKeyGen tests below without some seed premixing from
    // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
    return ROCKSDB_NAMESPACE::Hash64(key.data(), key.size(), raw_seed);
  }
  // For testing
  using KeyGen = StandardKeyGen;
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureAll;
  }
};

using TypesAndSettings_Coeff128 = DefaultTypesAndSettings;
struct TypesAndSettings_Coeff128Smash : public DefaultTypesAndSettings {
  static constexpr bool kUseSmash = true;
};
struct TypesAndSettings_Coeff64 : public DefaultTypesAndSettings {
  using CoeffRow = uint64_t;
};
struct TypesAndSettings_Coeff64Smash : public TypesAndSettings_Coeff64 {
  static constexpr bool kUseSmash = true;
};
struct TypesAndSettings_Coeff64Smash0 : public TypesAndSettings_Coeff64Smash {
  static constexpr bool kFirstCoeffAlwaysOne = false;
};

// Homogeneous Ribbon configurations
struct TypesAndSettings_Coeff128_Homog : public DefaultTypesAndSettings {
  static constexpr bool kHomogeneous = true;
  // Since our best construction success setting still has 1/1000 failure
  // rate, the best FP rate we test is 1/256
  using ResultRow = uint8_t;
  // Homogeneous only makes sense with sufficient slots for equivalent of
  // almost sure construction success
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnlyRare;
  }
};
struct TypesAndSettings_Coeff128Smash_Homog
    : public TypesAndSettings_Coeff128_Homog {
  // Smash (extra time to save space) + Homog (extra space to save time)
  // doesn't make much sense in practice, but we minimally test it
  static constexpr bool kUseSmash = true;
};
struct TypesAndSettings_Coeff64_Homog : public TypesAndSettings_Coeff128_Homog {
  using CoeffRow = uint64_t;
};
struct TypesAndSettings_Coeff64Smash_Homog
    : public TypesAndSettings_Coeff64_Homog {
  // Smash (extra time to save space) + Homog (extra space to save time)
  // doesn't make much sense in practice, but we minimally test it
  static constexpr bool kUseSmash = true;
};

// Less exhaustive mix of coverage, but still covering the most stressful case
// (only 50% construction success)
struct AbridgedTypesAndSettings : public DefaultTypesAndSettings {
  static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
    return kFailureOnly50Pct;
  }
};
struct TypesAndSettings_Result16 : public AbridgedTypesAndSettings {
  using ResultRow = uint16_t;
};
struct TypesAndSettings_Result32 : public AbridgedTypesAndSettings {
  using ResultRow = uint32_t;
};
struct TypesAndSettings_IndexSizeT : public AbridgedTypesAndSettings {
  using Index = size_t;
};
struct TypesAndSettings_Hash32 : public AbridgedTypesAndSettings {
  using Hash = uint32_t;
  static Hash HashFn(const Key& key, Hash raw_seed) {
    // This MurmurHash1 function does not pass tests below without the
    // seed premixing from StandardHasher. In fact, it needs more than
    // just a multiplication mixer on the ordinal seed.
    return ROCKSDB_NAMESPACE::Hash(key.data(), key.size(), raw_seed);
  }
};
struct TypesAndSettings_Hash32_Result16 : public AbridgedTypesAndSettings {
  using ResultRow = uint16_t;
};
struct TypesAndSettings_KeyString : public AbridgedTypesAndSettings {
  using Key = std::string;
};
struct TypesAndSettings_Seed8 : public AbridgedTypesAndSettings {
  // This is not a generally recommended configuration. With the configured
  // hash function, it would fail with SmallKeyGen due to insufficient
  // independence among the seeds.
  using Seed = uint8_t;
};
struct TypesAndSettings_NoAlwaysOne : public AbridgedTypesAndSettings {
  static constexpr bool kFirstCoeffAlwaysOne = false;
};
struct TypesAndSettings_AllowZeroStarts : public AbridgedTypesAndSettings {
  static constexpr bool kAllowZeroStarts = true;
};
struct TypesAndSettings_Seed64 : public AbridgedTypesAndSettings {
  using Seed = uint64_t;
};
//struct TypesAndSettings_Rehasher
//    : public StandardRehasherAdapter<AbridgedTypesAndSettings> {
//  using KeyGen = Hash64KeyGenWrapper<StandardKeyGen>;
//};
//struct TypesAndSettings_Rehasher_Result16 : public TypesAndSettings_Rehasher {
//  using ResultRow = uint16_t;
//};
//struct TypesAndSettings_Rehasher_Result32 : public TypesAndSettings_Rehasher {
//  using ResultRow = uint32_t;
//};
//struct TypesAndSettings_Rehasher_Seed64
//    : public StandardRehasherAdapter<TypesAndSettings_Seed64> {
//  using KeyGen = Hash64KeyGenWrapper<StandardKeyGen>;
//  // Note: 64-bit seed with Rehasher gives slightly better average reseeds
//};
//struct TypesAndSettings_Rehasher32
//    : public StandardRehasherAdapter<TypesAndSettings_Hash32> {
//  using KeyGen = Hash32KeyGenWrapper<StandardKeyGen>;
//};
//struct TypesAndSettings_Rehasher32_Coeff64
//    : public TypesAndSettings_Rehasher32 {
//  using CoeffRow = uint64_t;
//};
//struct TypesAndSettings_SmallKeyGen : public AbridgedTypesAndSettings {
//  // SmallKeyGen stresses the independence of different hash seeds
//  using KeyGen = SmallKeyGen;
//};
//struct TypesAndSettings_Hash32_SmallKeyGen : public TypesAndSettings_Hash32 {
//  // SmallKeyGen stresses the independence of different hash seeds
//  using KeyGen = SmallKeyGen;
//};
struct TypesAndSettings_Coeff32 : public DefaultTypesAndSettings {
  using CoeffRow = uint32_t;
};
struct TypesAndSettings_Coeff32Smash : public TypesAndSettings_Coeff32 {
  static constexpr bool kUseSmash = true;
};
struct TypesAndSettings_Coeff16 : public DefaultTypesAndSettings {
  using CoeffRow = uint16_t;
};
struct TypesAndSettings_Coeff16Smash : public TypesAndSettings_Coeff16 {
  static constexpr bool kUseSmash = true;
};

// For testing Poisson-distributed (or similar) statistics, get value for
// `stddevs_allowed` standard deviations above expected mean
// `expected_count`.
// (Poisson approximates Binomial only if probability of a trial being
// in the count is low.)
uint64_t PoissonUpperBound(double expected_count, double stddevs_allowed) {
  return static_cast<uint64_t>(
      expected_count + stddevs_allowed * std::sqrt(expected_count) + 1.0);
}

uint64_t PoissonLowerBound(double expected_count, double stddevs_allowed) {
  return static_cast<uint64_t>(std::max(
      0.0, expected_count - stddevs_allowed * std::sqrt(expected_count)));
}

uint64_t FrequentPoissonUpperBound(double expected_count) {
  // Allow up to 5.0 standard deviations for frequently checked statistics
  return PoissonUpperBound(expected_count, 5.0);
}

uint64_t FrequentPoissonLowerBound(double expected_count) {
  return PoissonLowerBound(expected_count, 5.0);
}

uint64_t InfrequentPoissonUpperBound(double expected_count) {
  // Allow up to 3 standard deviations for infrequently checked statistics
  return PoissonUpperBound(expected_count, 3.0);
}

uint64_t InfrequentPoissonLowerBound(double expected_count) {
  return PoissonLowerBound(expected_count, 3.0);
}

}  // namespace RIBBON_TEST

#endif  // ROCKSDB_RIBBON_TEST_H
