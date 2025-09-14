#pragma once

#include <iostream>
#include "common.h"

/// @brief 有理数
/// @details 分数をp/qの形で表現する
/// @note p,qは互いに素、q>0を保つ
struct rational
{
    long long p, q;
    void normalize()
    { // keep q positive
        if (q < 0)
            p *= -1, q *= -1;
        long long d = mathlib::gcd(p < 0 ? -p : p, q);
        if (d == 0)
            p = 0, q = 1;
        else
            p /= d, q /= d;
    }
    rational(long long p, long long q = 1) : p(p), q(q)
    {
        normalize();
    }
    rational &operator+=(const rational &a)
    {
        p = a.q * p + a.p * q;
        q = a.q * q;
        normalize();
        return *this;
    }
    rational &operator-=(const rational &a)
    {
        p = a.q * p - a.p * q;
        q = a.q * q;
        normalize();
        return *this;
    }
    rational &operator*=(const rational &a)
    {
        p *= a.p;
        q *= a.q;
        normalize();
        return *this;
    }
    rational &operator/=(const rational &a)
    {
        p *= a.q;
        q *= a.p;
        normalize();
        return *this;
    }
    rational &operator-()
    {
        p *= -1;
        return *this;
    }
    friend rational operator+(const rational &a, const rational &b) { return rational(a) += b; }
    friend rational operator*(const rational &a, const rational &b) { return rational(a) *= b; }
    friend rational operator-(const rational &a, const rational &b) { return rational(a) -= b; }
    friend rational operator/(const rational &a, const rational &b) { return rational(a) /= b; }
    friend bool operator<(const rational &a, const rational &b)
    { // avoid overflow
        return (long double)a.p * b.q < (long double)a.q * b.p;
    }
    friend bool operator<=(const rational &a, const rational &b) { return !(b < a); }
    friend bool operator>(const rational &a, const rational &b) { return b < a; }
    friend bool operator>=(const rational &a, const rational &b) { return !(a < b); }
    friend bool operator==(const rational &a, const rational &b) { return !(a < b) && !(b < a); }
    friend bool operator!=(const rational &a, const rational &b) { return (a < b) || (b < a); }
    friend std::ostream &operator<<(std::ostream &os, const rational &x)
    {
        printf("%.16f", (double)x.p / (double)x.q);
        return os;
    }
    friend std::istream &operator>>(std::istream &is, rational &x)
    {
        is >> x.p >> x.q;
        x.normalize();
        return is;
    }
};