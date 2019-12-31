#include<iostream>
typedef long long ll;
constexpr ll mod=1e9+7;
class mint {
 private:
  ll _num,_mod=mod;
  mint set(ll num){ 
      _num = num ;
      if(_num<0){
          if(_num>=-mod)_num=mod+_num;
          else _num=mod-(-_num)%mod;
      }
      else if(_num>=mod) _num%=mod;
      return *this;
  }
  ll imod()const{
    ll n=_mod-2;
    ll ans = 1,x=_num;
    while(n != 0){
        if(n&1) ans = ans*x%mod;
        x = x*x%mod;
        n = n >> 1;
    }
    return ans;
  }
 public:
  explicit mint(){ _num = 0; }
  explicit mint(ll num){
      _num = num;
      if(_num<0){
          if(_num>=-mod)_num=mod+_num;
          else _num=mod-(-_num)%mod;
      }
      else if(_num>=mod) _num%=mod;
  }
  explicit mint(ll num,ll M){
      _mod=M;
      _num=num;
      if(_num<0){
          if(_num>=-_mod)_num=mod+_num;
          else _num=mod-(-_num)%_mod;
      }
      else if(_num>=mod) _num%=_mod;
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
  mint operator/ (const ll x)const{ return *this/mint(x);}

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
struct rational
{
    long long p, q;
    long long gcd(ll x,ll y){return y?gcd(y,x%y):x;}
    void normalize()
    { // keep q positive
        if (q < 0)
            p *= -1, q *= -1;
        long long d = gcd(p < 0 ? -p : p, q);
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

template<typename T>
class Polynomial{//f(x)=a_0+a_1x+a_1x^2+...
    std::vector<T> fx;
  public:
    Polynomial(std::vector<T> fx={T(1)}):fx(fx){}
    T operator[](size_t k)const{return fx[k];}
    Polynomial operator=(const Polynomial &gx){
        fx.resize(gx.dim()+1);
        for(int i=0;i<fx.size();++i){
            fx[i]=gx[i];
        }
        return *this;
    }
    Polynomial operator+(const Polynomial &gx)const{
        size_t fs=std::max(fx.size(),gx.dim()+1);
        std::vector<T> hx(fs,T(0));
        for(int i=0;i<fs;++i){
            if(i<fx.size()) hx[i]+=fx[i];
            if(i<gx.dim()+1) hx[i]+=gx[i];
        }
        return hx;
    }
    Polynomial operator*(const Polynomial &gx)const{
        size_t fs=dim()+gx.dim()+1;
        std::vector<T> hx(fs,T(0));
        for(int i=0;i<fx.size();++i){
            for(int j=0;j<=gx.dim();++j){
                hx[i+j]+=fx[i]*gx[j];
            }
        }
        return hx;
    }
    Polynomial operator+=(const Polynomial &gx){*this=*this+gx;return *this;}
    Polynomial operator*=(const Polynomial &gx){*this=(*this)*gx;return *this;}
    Polynomial operator*(const T &x){
        for(int i=0;i<fx.size();++i){
            fx[i]*=x;
        }
        return *this;
    }
    size_t dim()const{return fx.size()-1;}
    T operator()(const T &x)const{
        T f=fx.back();
        for(int i=fx.size()-1;i>=0;--i){
            f=f*x+fx[i];
        }
        return f;
    }
    Polynomial integrate(){
        fx.push_back(0);
        for(int i=fx.size()-1;i>0;--i){
            fx[i]=fx[i-1]/i;
        }
        fx[0]=0;
        return *this;
    }
    Polynomial differencial(){
        for(int i=0;i+1<fx.size();++i){
            fx[i]=(i+1)*fx[i+1];
        }
        fx.pop_back();
        return *this;
    }
    T cintegrate(const T &a,const T &b){
        T d=fx.back()/dim(),u=fx.back()/dim();
        for(int i=fx.size()-1;i>0;--i){
            d=d*a+fx[i-1]/i;
            u=u*b+fx[i-1]/i;
        }
        return u-d;
    }
    friend std::ostream& operator<<(std::ostream &os, const Polynomial &gx){
        for(int i=0;i<gx.dim()+1;++i){
            os << gx[i] << "x^" <<i;
            if(i!=gx.dim()) os<<" + ";
        }
        return os;
    }
};
