#include<bits/stdc++.h>
#define REP(i,n) for(int i=0;i<(n);++i)
#define All(x) (x).begin(),(x).end()
using namespace std;
typedef long long ll;
typedef pair<ll,ll> P;
ll roundup(ll n,ll div){
    if(n>0) return (n-1)/div+1;
    else return n/div;    
}
P splitint(string n,int a){
    int Len=n.length();
    if(a<0||Len<a) return {-1,-1};
    string left,right;
    if(a!=0) left=n.substr(0,a);
    if(a!=Len) right=n.substr(a);
    return P(stoll(left),stoll(right));
}

bool kaibun(string s){//O(|S|)
    string  sdash=s;
    reverse(All(s));
    return s==sdash;
}
template<typename T>
class Ruiseki
{
  private:
    vector<T> LEFT, RIGHT;
    T d0;
    typedef function<T(T, T)> F;
    F f;
    F g;
    int N;

  public:
    Ruiseki(const vector<T> &a=vector<T>(), F f=[](T a,T b){return a+b;}, F g=[](T a,T b){rertun a-b;},T INI=0)
    :f(f),g(g),d0(INI)
    {
        N = a.size();
        LEFT.resize(N + 1);
        RIGHT.resize(N + 1);
        LEFT[0] = RIGHT[0] = INI;
        REP(i, N)
        {
            LEFT[i + 1] = f(LEFT[i], a[i]);
            RIGHT[i + 1] = f(RIGHT[i], a[N - i - 1]);
        }
    }
    T out(int l, int r)
    { //[l,r)の外の累積
        if(l>=r) return d0;
        if(l<0) return RIGHT[N-r];
        if(r>N) return LEFT[l];
        return F(LEFT[l], RIGHT[N - r]);
    }
    T in(int l,int r)
    { //[l,r)内の累積
        if(l>=r) return d0;
        return g(LEFT[r],LEFT[l]);
    }
};
