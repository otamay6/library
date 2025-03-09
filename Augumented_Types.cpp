#include<iostream>
typedef long long ll;
constexpr ll mod=1e9+7;
class mint {
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
    const double PI=acos(-1);
    std::vector<T> fx;
    using Complex=complex<double>;
    void DFT(std::vector<Complex> &F,int n,int sig=1)const{
        if(n==1) return;
        std::vector<Complex> f0(n/2),f1(n/2);
        for(int i=0;i<n/2;++i){
            f0[i]=F[2*i];
            f1[i]=F[2*i+1];
        }
        DFT(f0,n/2,sig);
        DFT(f1,n/2,sig);
        Complex z(cos(2.0*PI/n),sin(2.0*PI/n)*sig),zi=1;
        for(int i=0;i<n;++i){
            if(i<n/2) F[i]=f0[i]+zi*f1[i];
            else F[i]=f0[i-n/2]+zi*f1[i-n/2];
            zi*=z;
        }
        return;
    }
    void invDFT(std::vector<Complex> &f,int n)const{
        DFT(f,n,-1);
        for(int i=0;i<n;++i){
            f[i]/=n;
        }
        return;
    }
  public:
    Polynomial(const std::vector<T> &fx={T(1)}):fx(fx){}
    T operator[](size_t k)const{return fx[k];}
    size_t size()const{return fx.size();}
    size_t dim()const{return fx.size()-1;}
    T back()const{return fx.back();}
    Polynomial operator=(const Polynomial &gx){
        fx.resize(gx.size());
        for(int i=0;i<fx.size();++i){
            fx[i]=gx[i];
        }
        return *this;
    }

    Polynomial operator+(const Polynomial &gx)const{
        size_t fs=std::max(fx.size(),gx.size());
        std::vector<T> hx(fs,T(0));
        for(int i=0;i<fs;++i){
            if(i<fx.size()) hx[i]+=fx[i];
            if(i<gx.size()) hx[i]+=gx[i];
        }
        return Polynomial(hx);
    }
    Polynomial operator-(const Polynomial &gx)const{
        size_t fs=std::max(fx.size(),gx.size());
        std::vector<T> hx(fs,T(0));
        for(int i=0;i<fs;++i){
            if(i<fx.size()) hx[i]+=fx[i];
            if(i<gx.size()) hx[i]-=gx[i];
        }
        return Polynomial(hx);
    }
    
    Polynomial operator*(const Polynomial &gx)const{
        size_t fs=dim()+gx.size();
        int n=1;
        while(n<=fx.size()+gx.size()){
            n*=2;
        }
        std::vector<Complex> f(n,Complex(0)),g=f;
        int i_len=std::max(fx.size(),gx.size());
        for(int i=0;i<i_len;++i){
            if(i<fx.size()) f[i]=fx[i];
            if(i<gx.size()) g[i]=gx[i];
        }
        DFT(f,n);
        DFT(g,n);
        for(int i=0;i<n;++i) f[i]=f[i]*g[i];
        invDFT(f,n);
        std::vector<T> hx(fs);
        for(int i=0;i<fs;++i){
            if(T(f[i].real())!=f[i].real()) hx[i]=f[i].real()+0.5;
            else hx[i]=f[i].real();
        }
        return Polynomial(hx);
    }
    Polynomial operator/(const Polynomial &gx)const{
        std::vector<T> hx=fx;
        std::stack<T> st;
        while(hx.size()>=gx.size()){
            T t=hx.back()/gx.back();
            size_t n=hx.size()-gx.size();
            hx.pop_back();
            for(int i=0;i<gx.size()-1;++i){
                hx[i+n]-=t*gx[i];
            }
            st.push(t);
        }
        hx.clear();
        while(!st.empty()){
            hx.push_back(st.top());
            st.pop();
        }
        return Polynomial(hx);
    }
    Polynomial operator%(const Polynomial &gx)const{
        std::vector<T> hx=fx;
        while(hx.size()>=gx.size()){
            T t=hx.back()/gx.back();
            size_t n=hx.size()-gx.size();
            hx.pop_back();
            for(int i=0;i<gx.size()-1;++i){
                hx[i+n]-=t*gx[i];
            }
        }
        return Polynomial(hx);
    }
    Polynomial operator*(const T &x){
        for(int i=0;i<fx.size();++i){
            fx[i]*=x;
        }
        return *this;
    }
    Polynomial operator+=(const Polynomial &gx){*this=*this+gx;return *this;}
    Polynomial operator-=(const Polynomial &gx){*this=*this-gx;return *this;}
    Polynomial operator*=(const Polynomial &gx){*this=(*this)*gx;return *this;}
    Polynomial operator/=(const Polynomial &gx){*this=*this/gx;return *this;}
    Polynomial operator%=(const Polynomial &gx){*this=*this%gx;return *this;}

    T operator()(const T &x)const{
        if(fx.size()==0) return 0;
        if(fx.size()==1) return x; 
        T f=fx.back();
        for(int i=fx.size()-2;i>=0;--i){
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
    T cintegrate(const T &a,const T &b)const{
        if(fx.size()<=1) return 0;
        T d=fx.back()/size(),u=fx.back()/size();
        for(int i=fx.size()-1;i>0;--i){
            d=d*a+fx[i-1]/i;
            u=u*b+fx[i-1]/i;
        }
        d*=a;
        u*=b;
        return u-d;
    }
    friend std::ostream& operator<<(std::ostream &os, const Polynomial &gx){
        for(int i=0;i<gx.size();++i){
            os << gx[i] << "x^" <<i;
            if(i!=gx.dim()) os<<" + ";
        }
        return os;
    }
};
