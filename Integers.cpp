#include<vector>
#include<utility>
#include<algorithm>
#include<functional>
typedef long long ll;
ll npow(ll x, ll n){
    ll ans = 1;
    while(n != 0){
        if(n&1) ans = ans*x;
        x = x*x;
        n = n >> 1;
    }
    return ans;
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