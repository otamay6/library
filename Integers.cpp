#include<vector>
#include<utility>
#include<algorithm>
#include<functional>
typedef long long ll;
constexpr ll mod=1e9+7;
ll npow(ll x, ll n){
    ll ans = 1;
    if( x == 0 ) return 0;
    while(n != 0){
        if(n&1) ans = ans*x;
        x = x*x;
        n = n >> 1;
    }
    return ans;
}
ll mpow(ll x, ll n){ //x^n(mod) ←普通にpow(x,n)では溢れてしまうため，随時mod計算 2分累乗法だから早い
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
ll inv_mod(ll a){return mpow(a,mod-2);}
ll gcd(ll x,ll y){return y ? gcd(y,x%y) : x;};//xとyの最大公約数
ll lcm(ll x,ll y){return x/gcd(x,y)*y;}//xとyの最小公倍数
int digitsum(ll N,int a){
    if(N==0) return 0;
    int ret=0;
    ret+=digitsum(N/a,a)+N%a;
    return ret;
}
std::pair<ll,ll> bisearch(ll l,ll r,std::function<bool(ll)> f){
    while((l+1)<r){
        ll m=(l+r)/2;
        if(f(m)) r=m;
        else l=m;
    }
    return std::make_pair(l,r);
}
ll SQRT(ll n){if(n==1) return 1;return bisearch(0,n,[n](ll x){return x>n/x;}).first;}
ll manhattan(const std::pair<ll,ll> &a,const std::pair<ll,ll> &b){return llabs(a.first-b.first)+llabs(a.second-b.second);}
bool square(ll a){ll n=SQRT(a);return a==n*n;}
class Prime_Numbers{
 private:
    std::vector<int> PN;
    std::vector<bool> isp;
 public:
    Prime_Numbers(int N){
        isp.resize(N+1,true);
        isp[0]=isp[1]=false;
        for(int i=2;i<=N;++i) if(isp[i]){
            PN.push_back(i);
            for(int j=2*i;j<=N;j+=i) isp[j]=false;
        }
    }
    int operator[](int i){return PN[i];}
    int size(){return PN.size();}
    int back(){return PN.back();}
    bool isPrime(int q){return isp[q];}
};
class Divisor
{ //素因数分解をしてから約数列挙、分解結果はＰ(底,指数)でpfacにまとめている
private:
    std::vector<ll> F;
    std::vector<std::pair<int,int>> pfactorize;

public:
    Divisor(ll N){
        for(ll i=1; i*i<=N; i++){
            if(N%i == 0){
                F.push_back(i);
                if (i*i != N) F.push_back(N/i);
            }
        }
        sort(begin(F), end(F));
        Prime_Numbers p(SQRT(N) + 1);
        for(int i=0;i< p.size();++i){
            pfactorize.push_back(std::make_pair(p[i], 0));
            while(N%p[i] == 0){
                N /= p[i];
                pfactorize.back().second++;
            }
            if (pfactorize.size() == 0)
                pfactorize.pop_back();
        }
        if (N > 1)
            pfactorize.push_back(std::make_pair(N, 1));
    }
    int size() { return F.size(); }
    std::vector<std::pair<int,int>> pfac() { return pfactorize; }
    long long operator[](int k) { return F[k]; }
};
class Factorial
{ //階乗とその逆元を求めて計算に利用するクラス
  private:
    std::vector<ll> fac;
    std::vector<ll> ifac;
  public:
    Factorial(ll N)
    {
        fac.push_back(1);
        for(int i=0;i<N;++i)
        fac.push_back(fac[i] * (i + 1) % mod);
        ifac.resize(N + 1);
        ifac[N] = inv_mod(fac[N]);
        for(int i=0;i<N;++i)
        ifac[N - 1 - i] = (ifac[N - i] * (N - i)) % mod;
    }

    ll fact(ll a) { return fac[a]; }
    ll ifact(ll a) { return ifac[a]; }

    ll cmb(ll a, ll b)
    {
        if (a == 0 && b == 0)
            return 1;
        if (a < b || a < 0 || b < 0)
            return 0;
        ll tmp = ifact(a - b) * ifact(b) % mod;
        return tmp * fac[a] % mod;
    }
    ll per(ll a, ll b)
    {
        if (a == 0 && b == 0)
            return 1;
        if (a < b || a < 0 || b < 0)
            return 0;
        return fac[a] * ifac[a - b] % mod;
    }
};