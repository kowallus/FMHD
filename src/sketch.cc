#include "sketch.hpp"

#include <limits> /* max number which we can put into uint64_t */
#include "xxhash64.hpp" /* hashing function */
#include <unordered_map>
#include <cstring>
#include <iostream>
#include "commons.hpp"

// based on http://www.amsoftware.narod.ru/algo2.html
template <size_t K>
inline std::uint64_t maRushPrime1HashSimplified(const char *str) {
    std::uint64_t hash = K;
    std::uint32_t k;
    for (std::uint32_t j = 0; j < K/4; ) {
        memcpy(&k, str, 4);
        k += j++;
        hash ^= k;
        hash *= 171717;
        str += 4;
    }
    return (std::uint64_t)(hash);
}

template <size_t K>
inline std::uint64_t xxhash64(const char *str) {
    return XXHash64::hash(str, K, SEED);
}

template <size_t K, int HASH_ID>
inline uint64_t hash(char* input);

template <> inline uint64_t hash<20, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(char* input) { return maRushPrime1HashSimplified<20>(input); }
template <> inline uint64_t hash<21, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(char* input) { return maRushPrime1HashSimplified<20>(input); }
template <> inline uint64_t hash<24, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(char* input) { return maRushPrime1HashSimplified<24>(input); }
template <> inline uint64_t hash<20, XXHASH64_ID>(char* input) { return xxhash64<20>(input); }
template <> inline uint64_t hash<21, XXHASH64_ID>(char* input) { return xxhash64<21>(input); }
template <> inline uint64_t hash<24, XXHASH64_ID>(char* input) { return xxhash64<24>(input); }
template <> inline uint64_t hash<8, NONE_HASH_ID>(char* input) { return *((uint64_t*) input); }

template <size_t K, int HASH_ID>
inline void processInput(std::string &input, std::array<uint64_t, M> &hashes);

template <size_t K, int HASH_ID>
inline void processInputInBlocks(std::string &input, std::array<uint64_t, M> &hashes) {
    const size_t guard = input.size() - K + 1;
    int MULTI_K = 64;
    const size_t multi_guard = (guard / MULTI_K) * MULTI_K;
    uint64_t multiHashes[MULTI_K];
    size_t i = 0;
    for (; i < multi_guard; i += MULTI_K) {
        bool candidates = false;
        for (int j = 0; j < MULTI_K; j++) {
            multiHashes[j] = hash<K, HASH_ID>( input.data() + i + j );
//            if (~(multiHashes[j] & 0xFFC0000000000000))
            if ((multiHashes[j] >> 54) == 0)
                candidates = true;
        }
        if (candidates)
            for (int j = 0; j < MULTI_K; j++) {
                const uint64_t hashValue = multiHashes[j];
                uint64_t where = hashValue & MASK;
                hashes[where] = std::min(hashes[where], hashValue);
            }
    }
    for (; i < guard; ++i) {
        const uint64_t hashValue = hash<K, HASH_ID>( input.data() + i);
        uint64_t where = hashValue & MASK;
        hashes[where] = std::min(hashes[where], hashValue);
    }
}

template <size_t K, int HASH_ID>
inline void processInput(std::string &input, std::array<uint64_t, M> &hashes) {
    const size_t guard = input.size() - K + 1;
    for (size_t i = 0; i < guard; ++i) {
        const uint64_t hashValue = hash<K, HASH_ID>(input.data() + i);
        uint64_t where = hashValue & MASK;
        hashes[where] = std::min(hashes[where], hashValue);
    }
}

template<>
inline void processInput<8, NONE_HASH_ID>(std::string &input, std::array<uint64_t, M> &hashes)
{
    const size_t guard = input.size() - 4 + 1;
    for ( size_t i = 0; i<guard; ++i) {
        const uint64_t hashValue = *((uint64_t*) (input.data() + i));
        uint64_t where = hashValue & MASK;
        hashes[where] = std::min(hashes[where], hashValue);
    }
}

template<>
inline void processInputInBlocks<8, NONE_HASH_ID>(std::string &input, std::array<uint64_t, M> &hashes)
{
    const size_t guard = input.size() - 4 + 1;
    for ( size_t i = 0; i<guard; ++i) {
        const uint64_t hashValue = *((uint64_t*) (input.data() + i));
        uint64_t where = hashValue & MASK;
        hashes[where] = std::min(hashes[where], hashValue);
    }
}

template <size_t K, int HASH_ID>
inline std::array<uint64_t, M> get_sketch_full_template( std::string input )
{
    std::array<uint64_t, M> hashes{};
    hashes.fill(std::numeric_limits<uint64_t>::max());
    processInputInBlocks<K, HASH_ID>(input, hashes);
    reverseComplementInPlace(input);
    processInputInBlocks<K, HASH_ID>(input, hashes);
    return hashes;
}

template <int HASH_ID>
inline std::array<uint64_t, M> get_sketch_hash_template( std::string input, uint8_t kmerlen ) {
    if (kmerlen == 20) return get_sketch_full_template<20, HASH_ID>(input);
    if (kmerlen == 24) return get_sketch_full_template<24, HASH_ID>(input);
    if (HASH_ID != MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID && kmerlen == 21) return get_sketch_full_template<21, HASH_ID>(input);
    std::cerr << "unsupported kmer length: " << (int) kmerlen << " for hashID " << HASH_ID << std::endl;
    exit(EXIT_FAILURE);
}

std::array<uint64_t, M> get_sketch( std::string input, const uint8_t kmerlen, const int hash_id ) {
    if (hash_id == MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID)
        return get_sketch_hash_template<MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>(input, kmerlen);
    if (hash_id ==  XXHASH64_ID)
        return get_sketch_hash_template<XXHASH64_ID>(input, kmerlen);
    if (hash_id == NONE_HASH_ID)
        return get_sketch_full_template<8, NONE_HASH_ID>(input);
    std::cerr << "unsupported hash type id: " << hash_id << std::endl;
    exit(EXIT_FAILURE);
}

template std::array<uint64_t, M> get_sketch_hash_template<MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>( std::string , uint8_t );
template std::array<uint64_t, M> get_sketch_hash_template<XXHASH64_ID>( std::string , uint8_t );

template std::array<uint64_t, M> get_sketch_full_template<8, NONE_HASH_ID>( std::string );

template std::array<uint64_t, M> get_sketch_full_template<20, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>( std::string );
template std::array<uint64_t, M> get_sketch_full_template<21, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>( std::string );
template std::array<uint64_t, M> get_sketch_full_template<24, MA_RUSH_PRIME1_HASH_SIMPLIFIED_ID>( std::string );

template std::array<uint64_t, M> get_sketch_full_template<20, XXHASH64_ID>( std::string );
template std::array<uint64_t, M> get_sketch_full_template<21, XXHASH64_ID>( std::string );
template std::array<uint64_t, M> get_sketch_full_template<24, XXHASH64_ID>( std::string );
