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
    
    typedef union Double2Integer {
        double    Float;
        uint64_t  Integer;
    } Double2Integer;
    
    uint64_t ConvertDouble2Integer(double Decimal) {
        Double2Integer Integer = {.Float = Decimal};
        return Integer.Integer;
    }
    
#ifdef __cplusplus
}
#endif

#endif /* Ryu_Utilities_h */
