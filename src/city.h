/* reconstructed from VSEARCH authored by Torbjorn Rognes, Frederic Mahe and Tomas Flouri
 * the original source code can be found in https://github.com/torognes/vsearch
 * I have changed some behaviors and add some other functionality for better usage with this project
 */
#ifndef CITY_HASH_H_
#define CITY_HASH_H_

#include <stdlib.h>  // for size_t.
#include <stdint.h>
#include <string.h>
#include <utility>
#include <algorithm>

typedef std::pair<uint64_t, uint64_t> uint128;

inline uint64_t Uint128Low64(const uint128& x) { return x.first; }
inline uint64_t Uint128High64(const uint128& x) { return x.second; }

// Hash function for a byte array.
uint64_t CityHash64(const char *buf, size_t len);

// Hash function for a byte array.  For convenience, a 64-bit seed is also
// hashed into the result.
uint64_t CityHash64WithSeed(const char *buf, size_t len, uint64_t seed);

// Hash function for a byte array.  For convenience, two seeds are also
// hashed into the result.
uint64_t CityHash64WithSeeds(const char *buf, size_t len,
                           uint64_t seed0, uint64_t seed1);

// Hash function for a byte array.
uint128 CityHash128(const char *s, size_t len);

// Hash function for a byte array.  For convenience, a 128-bit seed is also
// hashed into the result.
uint128 CityHash128WithSeed(const char *s, size_t len, uint128 seed);

// Hash function for a byte array.  Most useful in 32-bit binaries.
uint32_t CityHash32(const char *buf, size_t len);

// Hash 128 input bits down to 64 bits of output.
// This is intended to be a reasonably good hash function.
inline uint64_t Hash128to64(const uint128& x) {
  // Murmur-inspired hashing.
  const uint64_t kMul = 0x9ddfea08eb382d69ULL;
  uint64_t a = (Uint128Low64(x) ^ Uint128High64(x)) * kMul;
  a ^= (a >> 47);
  uint64_t b = (Uint128High64(x) ^ a) * kMul;
  b ^= (b >> 47);
  b *= kMul;
  return b;
}

#endif  // CITY_HASH_H_
