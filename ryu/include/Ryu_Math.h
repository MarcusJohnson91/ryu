#include "Ryu_Tables.h"

#pragma once

#ifndef Ryu_Math_h
#define Ryu_Math_h

#ifdef __cplusplus
extern "C" {
#endif
    
    typedef enum RyuConstants {
        DOUBLE_POW5_INV_BITCOUNT = 122,
        DOUBLE_POW5_BITCOUNT     = 121,
        DOUBLE_MANTISSA_BITS     = 52,
        DOUBLE_EXPONENT_BITS     = 11,
        DOUBLE_BIAS              = 1023,
    } RyuConstants;
    
    uint64_t div5(const uint64_t x) {
        return x / 5;
    }
    
    uint64_t div10(const uint64_t x) {
        return x / 10;
    }
    
    uint64_t div100(const uint64_t x) {
        return x / 100;
    }
    
    uint64_t div1e8(const uint64_t x) {
        return x / 100000000;
    }
    
    // Returns e == 0 ? 1 : ceil(log_2(5^e)).
    int32_t pow5bits(const int32_t e) {
        // This approximation works up to the point that the multiplication overflows at e = 3529.
        // If the multiplication were done in 64 bits, it would fail at 5^4004 which is just greater
        // than 2^9297.
        assert(e >= 0);
        assert(e <= 3528);
        return (int32_t) (((((uint32_t) e) * 1217359) >> 19) + 1);
    }
    
    // Returns floor(log_10(2^e)).
    uint32_t log10Pow2(const int32_t e) {
        // The first value this approximation fails for is 2^1651 which is just greater than 10^297.
        assert(e >= 0);
        assert(e <= 1650);
        return (((uint32_t) e) * 78913) >> 18;
    }
    
    // Returns floor(log_10(5^e)).
    uint32_t log10Pow5(const int32_t e) {
        // The first value this approximation fails for is 5^2621 which is just greater than 10^1832.
        assert(e >= 0);
        assert(e <= 2620);
        return (((uint32_t) e) * 732923) >> 20;
    }
    
    uint64_t umul128(const uint64_t a, const uint64_t b, uint64_t *const productHi) {
        // The casts here help MSVC to avoid calls to the __allmul library function.
        const uint32_t aLo = (uint32_t)a;
        const uint32_t aHi = (uint32_t)(a >> 32);
        const uint32_t bLo = (uint32_t)b;
        const uint32_t bHi = (uint32_t)(b >> 32);
        
        const uint64_t b00 = (uint64_t)aLo  *bLo;
        const uint64_t b01 = (uint64_t)aLo  *bHi;
        const uint64_t b10 = (uint64_t)aHi  *bLo;
        const uint64_t b11 = (uint64_t)aHi  *bHi;
        
        const uint32_t b00Lo = (uint32_t)b00;
        const uint32_t b00Hi = (uint32_t)(b00 >> 32);
        
        const uint64_t mid1 = b10 + b00Hi;
        const uint32_t mid1Lo = (uint32_t)(mid1);
        const uint32_t mid1Hi = (uint32_t)(mid1 >> 32);
        
        const uint64_t mid2 = b01 + mid1Lo;
        const uint32_t mid2Lo = (uint32_t)(mid2);
        const uint32_t mid2Hi = (uint32_t)(mid2 >> 32);
        
        const uint64_t pHi = b11 + mid1Hi + mid2Hi;
        const uint64_t pLo = ((uint64_t)mid2Lo << 32) + b00Lo;
        
        *productHi = pHi;
        return pLo;
    }
    
    uint64_t shiftright128(const uint64_t lo, const uint64_t hi, const uint32_t dist) {
        // We don't need to handle the case dist >= 64 here (see above).
        assert(dist < 64);
        assert(dist > 0);
        return (hi << (64 - dist)) | (lo >> dist);
    }
    
    // Computes 5^-i in the form required by Ryu, and stores it in the given pointer.
    void double_computeInvPow5(const uint32_t i, uint64_t* const result) {
        const uint32_t base = (i + POW5_TABLE_SIZE - 1) / POW5_TABLE_SIZE;
        const uint32_t base2 = base * POW5_TABLE_SIZE;
        const uint32_t offset = base2 - i;
        const uint64_t* const mul = DOUBLE_POW5_INV_SPLIT2[base]; // 1/5^base2
        if (offset == 0) {
            result[0] = mul[0];
            result[1] = mul[1];
            return;
        }
        const uint64_t m = DOUBLE_POW5_TABLE[offset];
        uint64_t high1;
        const uint64_t low1 = umul128(m, mul[1], &high1);
        uint64_t high0;
        const uint64_t low0 = umul128(m, mul[0] - 1, &high0);
        const uint64_t sum = high0 + low1;
        if (sum < high0) {
            ++high1; // overflow into high1
        }
        // high1 | sum | low0
        const uint32_t delta = pow5bits(base2) - pow5bits(i);
        result[0] = shiftright128(low0, sum, delta) + 1 + ((POW5_INV_OFFSETS[i / 16] >> ((i % 16) << 1)) & 3);
        result[1] = shiftright128(sum, high1, delta);
    }
    
    // Computes 5^i in the form required by Ryu, and stores it in the given pointer.
    void double_computePow5(const uint32_t i, uint64_t *const result) {
        const uint32_t base = i / POW5_TABLE_SIZE;
        const uint32_t base2 = base * POW5_TABLE_SIZE;
        const uint32_t offset = i - base2;
        const uint64_t *const mul = DOUBLE_POW5_SPLIT2[base];
        if (offset == 0) {
            result[0] = mul[0];
            result[1] = mul[1];
            return;
        }
        const uint64_t m = DOUBLE_POW5_TABLE[offset];
        uint64_t high1;
        const uint64_t low1 = umul128(m, mul[1], &high1);
        uint64_t high0;
        const uint64_t low0 = umul128(m, mul[0], &high0);
        const uint64_t sum = high0 + low1;
        if (sum < high0) {
            ++high1; // overflow into high1
        }
        // high1 | sum | low0
        const uint32_t delta = pow5bits(i) - pow5bits(base2);
        result[0] = shiftright128(low0, sum, delta) + ((POW5_OFFSETS[base] >> offset) & 1);
        result[1] = shiftright128(sum, high1, delta);
    }
    
#ifdef __cplusplus
}
#endif

#endif /* Ryu_Math_h */
