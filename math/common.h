#pragma once

#include <functional>
#include <utility>
#include <cmath>

namespace mathlib {
    typedef long long ll;
    constexpr double eps=1e-9;
    const double PI=acos(-1);

    /// @brief nをdivで割ったときの切り上げ
    /// @param n 割られる数
    /// @param div 割る数
    /// @return n/divの切り上げ
    inline ll roundup(ll n,ll div){
        if(n>0) return (n-1)/div+1;
        else return n/div;    
    }

    inline ll npow(ll x, ll n){
        ll ans = 1;
        if( x == 0 ) return 0;
        while(n != 0){
            if(n&1) ans = ans*x;
            x = x*x;
            n = n >> 1;
        }
        return ans;
    }
    inline ll mpow(ll x, ll n, ll mod){ //x^n(mod) ←普通にpow(x,n)では溢れてしまうため，随時mod計算 2分累乗法だから早い
        ll ans = 1;
        x %= mod;
        if( x == 0 ) return 0;
        while(n != 0){
            if(n&1) ans = ans*x % mod;
            x = x*x % mod;
            n = n >> 1;
        }
        return ans;
    }
    inline ll inv_mod(ll a, ll mod){return mpow(a,mod-2, mod);}
    inline ll gcd(ll x,ll y){return y ? gcd(y,x%y) : x;};//xとyの最大公約数
    inline ll lcm(ll x,ll y){return x/gcd(x,y)*y;}//xとyの最小公倍数
    inline int digitsum(ll N,int a){
        if(N==0) return 0;
        int ret=0;
        ret+=digitsum(N/a,a)+N%a;
        return ret;
    }
    inline std::pair<ll,ll> bisearch(ll l,ll r,std::function<bool(ll)> f){
        while((l+1)<r){
            ll m=(l+r)/2;
            if(f(m)) r=m;
            else l=m;
        }
        return std::make_pair(l,r);
    }
    inline ll sqrt_floor(ll n){if(n==1) return 1;return bisearch(0,n,[n](ll x){return x>n/x;}).first;}
    inline ll manhattan_dist(const std::pair<ll,ll> &a,const std::pair<ll,ll> &b){return llabs(a.first-b.first)+llabs(a.second-b.second);}
    inline bool is_square(ll a){ll n=sqrt_floor(a);return a==n*n;}
}