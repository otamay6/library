#pragma once

#include <iostream>
#include "common.h"

template<mathlib::ll mod = 1e9 + 7>
class mint {
    using ll = mathlib::ll;
private:
  ll _num,_mod = mod;
  mint set(ll num){ 
      _num = num ;
      if(_num<0){
          if(_num>=-_mod)_num=_mod+_num;
          else _num=_mod-(-_num)%_mod;
      }
      else if(_num>=_mod) _num%=_mod;
      return *this;
  }
  ll imod()const{
    ll n=_mod-2;
    ll ans = 1,x=_num;
    while(n != 0){
        if(n&1) ans = ans*x%_mod;
        x = x*x%_mod;
        n = n >> 1;
    }
    return ans;
  }
 public:
  explicit mint(){ _num = 0; }
  explicit mint(ll num){
      _num = num;
      if(_num<0){
          if(_num>=-_mod)_num=_mod+_num;
          else _num=_mod-(-_num)%_mod;
      }
      else if(_num>=_mod) _num%=_mod;
  }
  explicit mint(ll num,ll M){
      _mod=M;
      _num=num;
      if(_num<0){
          if(_num>=-_mod)_num=_mod+_num;
          else _num=_mod-(-_num)%_mod;
      }
      else if(_num>=_mod) _num%=_mod;
  }
  mint(const mint &cp){_num=cp._num;_mod=cp._mod;}
  
  mint operator+ (const mint &x)const{ return mint(_num + x._num , _mod); }
  mint operator- (const mint &x)const{ return mint(_num - x._num , _mod);}
  mint operator* (const mint &x)const{ return mint(_num * x._num , _mod); }
  mint operator/ (const mint &x)const{ return mint(_num * x.imod() , _mod);}
  
  mint operator+=(const mint &x){ return set(_num + x._num); }
  mint operator-=(const mint &x){ return set(_num - x._num); }
  mint operator*=(const mint &x){ return set(_num * x._num); }
  mint operator/=(const mint &x){ return set(_num * x.imod());}

  mint operator= (const ll x){ return set(x); }
  mint operator+ (const ll x)const{return *this + mint(x,_mod); }
  mint operator- (const ll x)const{ return *this - mint(x,_mod); }
  mint operator* (const ll x)const{ return *this * mint(x,_mod); }
  mint operator/ (const ll x)const{ return *this/mint(x, _mod);}

  mint operator+=(const ll x){ *this = *this + x;return *this; }
  mint operator-=(const ll x){ *this = *this - x;return *this; }
  mint operator*=(const ll x){ *this = *this * x;return *this;}
  mint operator/=(const ll x){ *this = *this / x;return *this;}

  bool operator==(const mint &x)const{return _num==x._num;}
  bool operator!=(const mint &x)const{return _num!=x._num;}

  friend mint operator+(ll x,const mint &m){return mint(m._num + x , m._mod);}
  friend mint operator-(ll x,const mint &m){return mint( x - m._num , m._mod);}
  friend mint operator*(ll x,const mint &m){return mint(m._num * (x % m._mod) , m._mod);}
  friend mint operator/(ll x,const mint &m){return mint(m.imod() * (x % m._mod) , m._mod);}

  explicit operator ll() { return _num; }
  explicit operator int() { return (int)_num; }
  
  friend std::ostream& operator<<(std::ostream &os, const mint &x){ os << x._num; return os; }
  friend std::istream& operator>>(std::istream &is, mint &x){ll val; is>>val; x.set(val); return is;}
};