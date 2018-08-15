/*  Written in 2017 by Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>.

   NOTE: as of 2017-10-08, this generator has a different multiplier (a
   fixed-point representation of the golden ratio), which eliminates
   linear dependencies from one of the lowest bits. The previous
   multiplier was 1181783497276652981 (M_8 in the paper). If you need to
   tell apart the two generators, you can refer to this generator as
   xorshift1024φ and to the previous one as xorshift1024*M_8.

   This is a fast, high-quality generator. If 1024 bits of state are too
   much, try a xoroshiro128+ generator.

   Note that the two lowest bits of this generator are LFSRs of degree
   1024, and thus will fail binary rank tests. The other bits needs a much
   higher degree to be represented as LFSRs.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s.

   See URL: https://github.com/jj1bdx/xorshiftplus-c/blob/original-20160106/xorshift1024star.c
*/

#include <stdint.h>
#include <string.h>



uint64_t xrs[16];
int xrp = 0;


uint64_t xorShift1024(void)
{
  const uint64_t s0 = xrs[xrp];
  uint64_t s1 = xrs[xrp = (xrp + 1) & 15];
  s1 ^= s1 << 31;
  xrs[xrp] = s1 ^ s0 ^ (s1 >> 11) ^ (s0 >> 30);
  return xrs[xrp] * UINT64_C(1181783497276652981);
}

