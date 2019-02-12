#include "Ryu_Math.h"

#pragma once

#ifndef Ryu_Utilities_h
#define Ryu_Utilities_h

#ifdef __cplusplus
extern "C" {
#endif
    
    int copy_special_str(char *const result, const bool sign, const bool exponent, const bool mantissa) {
        if (mantissa) {
            memcpy(result, "NaN", 3);
            return 3;
        }
        if (sign) {
            result[0] = '-';
        }
        if (exponent) {
            memcpy(result + sign, "Infinity", 8);
            return sign + 8;
        }
        memcpy(result + sign, "0E0", 3);
        return sign + 3;
    }
    
    uint64_t double_to_bits(const double d) {
        uint64_t bits = 0;
        memcpy(&bits, &d, sizeof(double));
        return bits;
    }
    
#ifdef __cplusplus
}
#endif

#endif /* Ryu_Utilities_h */
