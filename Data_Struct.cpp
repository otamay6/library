#include<vector>
#include<utility>
#include<functional>
template<typename T> class SegmentTree{
  private:
    typedef std::function<T(T,T)> F;
    int n;
    T d0;
    std::vector<T> vertex;
    F f;
    F g;
  public:
    SegmentTree(F f,F g,T d):d0(d),f(f),g(g){}
    void init(int _n){
        n=1;
        while(n<_n) n*=2;
        vertex.resize(2*n-1,d0);
    }
    void build(const std::vector<T> &v){
        int n_=v.size();
        init(n_);
        for(int i=0;i<n_;i++) vertex[n+i-1]=v[i];
        for(int i=n-2;i>=0;i--)
        vertex[i]=f(vertex[2*i+1],vertex[2*i+2]);
    }
    void update(int i,T x){
        int k=i+n-1;
        vertex[k]=g(vertex[k],x);
        while(k>0){
            k=(k-1)/2;
            vertex[k]=f(vertex[2*k+1],vertex[2*k+2]);
        }
        return;
    }
    T query(int l,int r){
        T vl=d0,vr=d0;
        l+=n-1;
        r+=n-1;
        for(;l<=r;l/=2,r=r/2-1){
            if(l%2==0) vl=f(vl,vertex[l]);
            if(r&1) vr=f(vr,vertex[r]);
        }
        return f(vl,vr);
    }
};
class UnionFind
{ //UnionFind木
private:
    std::vector<int> Parent, es;
    std::vector<long long> diff_weight;

public:
    UnionFind(int N)
    {
        es.resize(N, 0);
        Parent.resize(N, -1);
        diff_weight.resize(N, 0LL);
    }

    int root(int A)
    {
        if (Parent[A] < 0)
            return A;
        else
        {
            int r = root(Parent[A]);
            diff_weight[A] += diff_weight[Parent[A]]; // 累積和をとる
            return Parent[A] = r;
        }
    }
    bool issame(int A, int B)
    {
        return root(A) == root(B);
    }
    long long weight(int x)
    {
        root(x); // 経路圧縮
        return diff_weight[x];
    }
    long long diff(int x, int y)
    {
        return weight(y) - weight(x);
    }
    int size(int A)
    {
        return -Parent[root(A)];
    }
    int eize(int A)
    {
        return es[root(A)];
    }
    bool connect(int A, int B)
    {
        A = root(A);
        B = root(B);
        if (A == B)
            return false;
        if (size(A) < size(B))
            std::swap(A, B);
        Parent[A] += Parent[B];
        es[A] += es[B] + 1;
        Parent[B] = A;
        return true;
    }
    void unite(int A, int B)
    {
        A = root(A);
        B = root(B);
        if (A == B)
        {
            es[A]++;
            return;
        }
        if (size(A) < size(B))
            std::swap(A, B);
        Parent[A] += Parent[B];
        es[A] += es[B] + 1;
        Parent[B] = A;
        return;
    }
    bool merge(int A, int B, long long w)
    {
        // x と y それぞれについて、 root との重み差分を補正
        w += weight(A);
        w -= weight(B);
        A = root(A);
        B = root(B);
        if (A == B)
            return false;
        if (size(A) < size(B))
            std::swap(A, B), w = -w;
        Parent[A] += Parent[B];
        Parent[B] = A;
        // x が y の親になるので、x と y の差分を diff_weight[y] に記録
        diff_weight[B] = w;
        return true;
    }
};

template<typename T>
class DAG{
  private:
    int v;
    std::vector<std::vector<std::pair<int,T>>> to;
    std::vector<int> rank;
    std::vector<T> dp;
  public:
    DAG(int v,const T &x):v(v){
        to.resize(v);
        rank.resize(v);
        dp.resize(v,x);
        std::iota(rank.begin(),rank.end(),0);
    }
    void add(int a,int b,const T &c=0){
        to[a].push_back(std::make_pair(b,c));
        if(rank[a]>rank[b]) std::swap(rank[a],rank[b]);
    }
    void set(int i,const T &x){
        dp[i]=x;
    }
    std::vector<T> calc(std::function<T(T,T)> f){
        std::vector<int> vertex(v);
        std::vector<T> tmp=dp;
        for(int i=0;i<v;++i) vertex[rank[i]]=i;
        for(int i=0;i<v;++i){
            for(auto e:to[vertex[i]]){
                tmp[e.first]=f(tmp[e.first],tmp[i]+e.second);
            }
        }
        return tmp;
    }
};

template<typename T>
class FFT{
  private:
    using Complex = std::complex<double>;
    std::vector<Complex> C;
    void DFT(std::vector<Complex> &F,int n,int sig=1){
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
    void invDFT(std::vector<Complex> &f,int n){
        DFT(f,n,-1);
        for(int i=0;i<n;++i){
            f[i]/=n;
        }
        return;
    }
  public:
    FFT(const std::vector<T> &A,const std::vector<T> &B){
        int n=1;
        while(n<=A.size()+B.size()){
            n*=2;
        }
        std::vector<Complex> h(n,Complex(0)),g(n,Complex(0));
        for(int i=0;i<n;++i){
            if(i<A.size()) h[i]=Complex(A[i]);
            if(i<B.size()) g[i]=Complex(B[i]);
        }
        DFT(h,n);
        DFT(g,n);
        C.resize(n);
        for(int i=0;i<n;++i) C[i]=h[i]*g[i];
        invDFT(C,n);
        for(int i=0;i<n;++i) if(T(C[i].real())!=C[i].real()){
            C[i]=Complex(C[i].real()+0.5);
        }
    }
    T operator[](int k)const{return C[k].real();}
};