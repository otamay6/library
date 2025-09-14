#pragma once
#include "common.h"

template<mathlib::ll mod = 1e9 + 7>
class ModuloFactorial
{ //階乗とその逆元を求めて計算に利用するクラス
  private:
    std::vector<mathlib::ll> fac;
    std::vector<mathlib::ll> ifac;
  public:
    ModuloFactorial(mathlib::ll N)
    {
        fac.push_back(1);
        for(int i=0;i<N;++i)
            fac.push_back(fac[i] * (i + 1) % mod);
        ifac.resize(N + 1);
        ifac[N] = mathlib::inv_mod(fac[N], mod);
        for(int i=0;i<N;++i)
            ifac[N - 1 - i] = (ifac[N - i] * (N - i)) % mod;
    }

    mathlib::ll fact(mathlib::ll a) { return fac[a]; }
    mathlib::ll ifact(mathlib::ll a) { return ifac[a]; }

    mathlib::ll cmbination(mathlib::ll a, mathlib::ll b)
    {
        if (a == 0 && b == 0)
            return 1;
        if (a < b || a < 0 || b < 0)
            return 0;
        mathlib::ll tmp = ifact(a - b) * ifact(b) % mod;
        return tmp * fac[a] % mod;
    }
    mathlib::ll per(mathlib::ll a, mathlib::ll b)
    {
        if (a == 0 && b == 0)
            return 1;
        if (a < b || a < 0 || b < 0)
            return 0;
        return fac[a] * ifac[a - b] % mod;
    }
};