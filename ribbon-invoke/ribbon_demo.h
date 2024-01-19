//
// Created by Rain Night on 2023/12/28.
//

#ifndef ROCKSDB_RIBBON_COMPILE_BLOOM_H
#define ROCKSDB_RIBBON_COMPILE_BLOOM_H

#include <stdint.h>
#include "ribbon-lib/ribbon_config.h"
#include "ribbon-lib/ribbon_impl.h"
#include "ribbon-lib/slice.h"
#include "ribbon-lib/hash.h"

namespace RIBBON_FILTER {

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
        using Key = Slice;
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
        static const std::vector<ConstructionFailureChance>& FailureChanceToTest() {
            return kFailureAll;
        }
    };

    struct Settings_Coeff128_Homog : public DefaultTypesAndSettings {
        static constexpr bool kHomogeneous = true;
        static constexpr bool kUseSmash = false;

        static const std::vector<ConstructionFailureChance> &FailureChanceToTest() {
            return kFailureOnlyRare;
        }
    };

    struct tmp_Settings : public DefaultTypesAndSettings {
        static constexpr bool kHomogeneous = true;
        static constexpr bool kUseSmash = true;

        static const std::vector<ConstructionFailureChance> &FailureChanceToTest() {
            return kFailureOnlyRare;
        }
    };

    struct tmp_Settings2 : public DefaultTypesAndSettings{
        using Key = Hash;
        static constexpr bool kHomogeneous = true;
        static constexpr bool kUseSmash = true;
        static Hash HashFn(const Hash& key, uint64_t raw_seed) {
            return (key ^ raw_seed) * kRehashFactor;
        }

        static Hash HashFn2(const Slice& key, uint64_t raw_seed) {
            // This version 0.7.2 preview of XXH3 (a.k.a. XXPH3) function does
            // not pass SmallKeyGen tests below without some seed premixing from
            // StandardHasher. See https://github.com/Cyan4973/xxHash/issues/469
            return ROCKSDB_NAMESPACE::Hash64(key.data(), key.size(), raw_seed);
        }
    private:
        static constexpr Hash kRehashFactor =
                static_cast<Hash>(0x6193d459236a3a0dULL);
    };

} // namespace
#endif