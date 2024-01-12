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


}  // namespace RIBBON_TEST

#endif  // ROCKSDB_RIBBON_TEST_H
