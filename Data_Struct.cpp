#include<vector>
#include<queue>
#include<stack>
#include<utility>
#include<functional>
#include<random>

/// @brief セグ木
/// @tparam T モノイド
template<typename T> class SegmentTree{
  private:
    typedef std::function<T(T,T)> F;
    int n;
    T d0;
    std::vector<T> vertex;
    F f;
    F g;
  public:
    /// @brief セグ木作成
    /// @param f queryの適用関数
    /// @param g updateの適用関数
    /// @param d 零元
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
    /// @brief [l,r]のqueryを求める
    /// @param l 
    /// @param r 
    /// @return 
    T query(int l,int r){
        T vl=d0,vr=d0;
        l+=n-1;
        r+=n-1;
        for(;l<=r;l/=2,r=r/2-1){
            if(~l&1) vl=f(vl,vertex[l]);
            if(r&1) vr=f(vertex[r],vr);
        }
        return f(vl,vr);
    }
};


template<typename T>
class dynamic_segtree{
  private:
    using F = std::function<T(T,T)>;
    F f;
    F g;
    T ei;
    struct node_t{
        int L,R;
        T val;
        node_t *ch[2];
        node_t(int L,int R,T val):L(L),R(R),val(val){
            ch[0]=NULL;
            ch[1]=NULL;
        }
    };
    node_t *root;
    node_t *update(node_t* nd,int k,T x){
        if(nd->L==nd->R){
            nd->val = g(nd->val,x);
            return nd;
        }
        if(k<=(nd->R+nd->L)/2){
            if(!nd->ch[0]) nd->ch[0] = new node_t(nd->L,(nd->L+nd->R)/2,ei);
            nd->ch[0] = update(nd->ch[0],k,x);
        }
        else{
            if(!nd->ch[1]) nd->ch[1] = new node_t((nd->L+nd->R)/2+1,nd->R,ei);
            nd->ch[1] = update(nd->ch[1],k,x);
        }
        nd->val = f(nd->ch[0]?nd->ch[0]->val:ei,nd->ch[1]?nd->ch[1]->val:ei);
        return nd;
    }
    T query(node_t *nd,int l,int r){
        if(nd->R<l||r<nd->L) return ei;
        if(l<=nd->L&&nd->R<=r) return nd->val;
        return f(nd->ch[0]?query(nd->ch[0],l,r):ei,nd->ch[1]?query(nd->ch[1],l,r):ei);
    }
  public:
    dynamic_segtree(F f,T ei):f(f),ei(ei){}
    ~dynamic_segtree(){delete root;root = nullptr;}
    void init(int n){
        root = new node_t(0,n-1,ei);
    }
    void update(int k,T x,F fc){
        g=fc;
        update(root,k,x);
    }
    T query(int l,int r){
        return query(root,l,r);
    }
};

template<typename T,typename E>
class LazySegmentTree{
    private:
        using F = std::function<T(T,T)>;
        using G = std::function<T(T,E)>;
        using H = std::function<E(E,E)>;
        using HH = std::function<E(E)>;
        int n;
        std::vector<T> dat;
        std::vector<E> laz;
        T t0;
        E e0;
        F f;
        G g;
        H h;
        HH hh;
        void eval(int k){
            if(laz[k]==e0) return;
            if(k<n-1){
                laz[2*k+1]=h(laz[2*k+1],hh(laz[k]));
                laz[2*k+2]=h(laz[2*k+2],hh(laz[k]));
            }
            dat[k]=g(dat[k],laz[k]);
            laz[k]=e0;
        }
        void merge_dat(int k){
            eval(2*k+1);
            eval(2*k+2);
            dat[k]=f(dat[2*k+1],dat[2*k+2]);
        }
    public:
        LazySegmentTree(F f,G g,H h,T t0,E e0,HH hh=[](E x){return x;}):f(f),g(g),h(h),t0(t0),e0(e0),hh(hh){}
        void init(int n_){
            n=1;
            while(n<n_){
                n*=2;
            }
            dat.resize(2*n-1,t0);
            laz.resize(2*n-1,e0);
        }
        void build(const std::vector<T> &v){
            int n_=v.size();
            init(n_);
            for(int i=0;i<n_;i++) dat[n+i-1]=v[i];
            for(int i=n-2;i>=0;i--)
            dat[i]=f(dat[2*i+1],dat[2*i+2]);
        }
        void _update(int a,int b,E x,int k,int l,int r){
            if(b<l||r<a||k>=2*n-1) return;
            eval(k);
            if(a<=l&&r<=b){
                laz[k]=h(laz[k],x);
                eval(k);
            }
            else{
                _update(a,b,x,2*k+1,l,(l+r)/2);
                _update(a,b,x,2*k+2,(l+r)/2+1,r);
                dat[k]=f(dat[2*k+1],dat[2*k+2]);
            }
        }
        void update(int l,int r,E x){
            l+=n-1,r+=n-1;
            for(;l<=r;l>>=1,r=r/2-1){
                if(~l&1){
                    laz[l]=h(laz[l],x);
                    int k=l;
                    while(k>0&&~k&1){
                        eval(k);
                        k=k/2-1;
                        merge_dat(k);
                    }
                    eval(k);
                }
                if(r&1){
                    laz[r]=h(laz[r],x);
                    int k=r;
                    while(k&1){
                        eval(k);
                        k>>=1;
                        merge_dat(k);
                    }
                    eval(k);
                }
                x=h(x,x);
            }
            while(l>0){
                l=(l-1)/2;
                merge_dat(l);
            }
        }
        T _query(int a,int b,int k,int l,int r){
            if(b<l||r<a||k>=2*n-1) return t0;
            eval(k);
            if(a<=l&&r<=b) return dat[k];
            return f(_query(a,b,2*k+1,l,(l+r)/2),_query(a,b,2*k+2,(l+r)/2+1,r));
        }
        T query(int l,int r){
            return _query(l,r,0,0,n-1);
        }
        void set_val(int i,T x){
        int k=i+n-1;
        dat[k]=x;
        laz[k]=e0;
        while(k>0){
            k=(k-1)/2;
            dat[k]=f(dat[2*k+1],dat[2*k+2]);
        }
        return;
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

class HLD{
  private:
    using Graph=std::vector<std::vector<int>>;
    int root,N;
    std::vector<int> size,par,Restore,newid,front,next,Depth,leaf;
    Graph graph;
    int dfs(int u){
        int cnt=1,MAX=0;
        for(auto v:graph[u]){
            if(v==par[u]) continue;
            par[v]=u;
            cnt+=dfs(v);
            if(size[v]>MAX){
                MAX=size[v];
                next[u]=v;
            }
        }
        return size[u]=cnt;
    }
    void hld(int r){
        std::stack<int> que;
        int id=0;
        que.push(r);
        while(!que.empty()){
            int u=que.top();que.pop();
            newid[u]=id;
            Restore[id++]=u;
            front[u]=u;
            while(next[u]!=-1){
                for(auto v:graph[u]){
                    if(v==par[u]) continue;
                    if(v!=next[u]){
                        Depth[v]=Depth[u]+1;
                        que.push(v);
                    }
                    else{
                        newid[v]=id;
                        Restore[id++]=v;
                        front[v]=front[u];
                        Depth[v]=Depth[u];
                    }
                }
                u=next[u];
            }
        }
    }
  public:
    HLD(Graph &g,int root=0):graph(g),root(root),N(g.size()),
    size(N),par(N,-1),Restore(N),newid(N),front(N),next(N,-1),Depth(N),leaf(N,-1){
        dfs(root);
        hld(root);
    }
    int lca(int u,int v)const{
        while(front[u]!=front[v]){
            if(Depth[u]>Depth[v]) std::swap(u,v);
            v = par[front[v]];
        }
        return newid[u]<newid[v]?u:v;
    }
    std::vector<std::pair<int,int>> vquery(int u,int v)const{
        std::vector<std::pair<int,int>> res;
        int rt=lca(u,v);
        while(Depth[rt]<Depth[u]){
            res.emplace_back(newid[front[u]],newid[u]);
            u = par[front[u]];
        }
        while(Depth[rt]<Depth[v]){
            res.emplace_back(newid[front[v]],newid[v]);
            v = par[front[v]];
        }
        res.emplace_back(newid[rt],std::max(newid[u],newid[v]));
        return res;
    }
    std::vector<std::pair<int,int>> equery(int u,int v)const{
        //頂点idから親に向かう辺の番号をid-1とする
        std::vector<std::pair<int,int>> res;
        int rt=lca(u,v);
        while(Depth[rt]<Depth[u]){
            res.emplace_back(newid[front[u]]-1,newid[u]-1);
            u = par[front[u]];
        }
        while(Depth[rt]<Depth[v]){
            res.emplace_back(newid[front[v]]-1,newid[v]-1);
            v = par[front[v]];
        }
        int R = std::max(newid[u],newid[v]);
        if(newid[rt]!=R) res.emplace_back(newid[rt],R-1);
        return res;
    }
    std::pair<int,int> tquery(int u){
        if(leaf[u]!=-1) return std::make_pair(newid[u],leaf[u]);
        leaf[u]=newid[u];
        for(auto v:graph[u]){
            if(v==par[u]) continue;
            tquery(v);
            leaf[u]=std::max(leaf[u],leaf[v]);
        }
        return std::make_pair(newid[u],leaf[u]);
    }
    int id(int u)const{return newid[u];}
    int restore(int ID)const{return Restore[ID];}
    int depth(int u)const{return Depth[u];}
    int parent(int u)const{return par[u];}
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
    void DFT(std::vector<Complex> &F,size_t n,int sig=1){
        if(n==1) return;
        std::vector<Complex> f0(n/2),f1(n/2);
        for(size_t i=0;i<n/2;++i){
            f0[i]=F[2*i];
            f1[i]=F[2*i+1];
        }
        DFT(f0,n/2,sig);
        DFT(f1,n/2,sig);
        Complex z(cos(2.0*PI/n),sin(2.0*PI/n)*sig),zi=1;
        for(size_t i=0;i<n;++i){
            if(i<n/2) F[i]=f0[i]+zi*f1[i];
            else F[i]=f0[i-n/2]+zi*f1[i-n/2];
            zi*=z;
        }
        return;
    }
    void invDFT(std::vector<Complex> &f,size_t n){
        DFT(f,n,-1);
        for(size_t i=0;i<n;++i){
            f[i]/=n;
        }
        return;
    }
  public:
    FFT(const std::vector<T> &A,const std::vector<T> &B){
        size_t n=1;
        while(n<=A.size()+B.size()){
            n*=2;
        }
        std::vector<Complex> g(n,Complex(0));
        C.resize(n,Complex(0));
        size_t i_len=std::max(A.size(),B.size());
        for(size_t i=0;i<i_len;++i){
            if(i<A.size()) C[i]=Complex(A[i]);
            if(i<B.size()) g[i]=Complex(B[i]);
        }
        DFT(C,n);
        DFT(g,n);
        for(size_t i=0;i<n;++i) C[i]=C[i]*g[i];
        invDFT(C,n);
        for(size_t i=0;i<n;++i) if(T(C[i].real())!=C[i].real()){
            C[i]=Complex(C[i].real()+0.5);
        }
    }
    T operator[](int k)const{return T(C[k].real());}
};

/// @brief Treap
/// @tparam T : キーの型
template<class T>
class Treap {
    std::random_device rd;
    std::mt19937 mt;

    struct Node{
        T key;
        int priority;
        Node *l;
        Node *r;
        Node(T key, int priority)
        :key(key),priority(priority), l(nullptr),r(nullptr) {

        }
    };
    Node *root;
    void deconstruct(Node *root) {
        if(root->l) deconstruct(root->l);
        if(root->r) deconstruct(root->r);
        delete root;
    }
    /// @brief Nodeをkeyで左右に分割する
    /// @param root Treapの根
    /// @param key 
    /// @param l 左の部分木
    /// @param r 右の部分木
    void split(Node *root, T key, Node **l, Node **r) {
        if(!root) {
            *l = *r = nullptr;
        }
        else if(key < root->key) {
            split(root->l, key, l, root->l);
            *r = root;
        }
        else {
            split(root->r, key, root->r, r);
            *l = root;
        }
    }

    void insert(Node **root, Node *item) {
        if(!root) {
            *root = item;
        }
        else if(item->priority > root->priority){
            split(root, item->key, &(item->l), &(item->r));
            *root = item;
        } else {
            if(item->key < root->key) {
                insert(root->l, item);
            }
            else {
                insert(root->r, item);
            }
        }
    }

    void merge(Node **root, Node *left, Node *right) {
        if(!left) {
            *root = right;
        }
        else if(!right) {
            *root = left;
        }
        else if(left->priority, right->priority) {
            merge(left->r, left->r, right);
            *root = left;
        }
        else {
            merge(right->l, left, right->l);
            *root = right;
        }
    }

    void erase(Node **root, T key) {
        if(*root->key == key) {
            merge(root, *root->l, *root->r);
        }
        else {
            if(key < *root->key) {
                erase(*root->l, key);
            }
            else {
                erase(*root->r, key);
            }
        }
    }

    bool find(Node *root, T key) {
        if(!root) {
            return false;
        }
        else if(root->key == key) {
            return true;
        }
        else {
            if(key < root->key){
                return find(root->l, key);
            }
            else {
                return find(root->r, key);
            }
        }
    }

public:
    Treap() :mt(rd()), root(nullptr){};
    ~Treap() { deconstruct(root); }
    /// @brief keyを持つノードを挿入する
    /// @param key 
    void insert(T key) {
        insert(&root, new Node(key, mt()));
    }

    /// @brief keyを持つノードを1つ削除する
    /// @param key 
    void erase(T key) {
        erase(&root, key);
    }

    /// @brief keyが存在するか判定する
    /// @param key 
    /// @return 
    bool find(T key) {
        return find(root, key);
    }
};

/// @brief 操作をO(logN)でできる配列(insert, erase, reverse, rotate, range min, range add)
/// @tparam T : 配列要素の型
template<class T>
class ImplicitTreap {
    std::random_device rd;
    std::mt19937 mt;
    T INF;
    struct Node {
        T value, min, lazy;
        int priority, cnt;
        bool rev;
        Node *l, *r;
        Node(int value, int priority, T INF)
        : value(value), min(INF), lazy(0), priority(priority), cnt(1), rev(false), l(nullptr), r(nullptr) 
        {

        }
    } *root = nullptr;
    using Tree = Node *;

    int cnt(Tree t) {
        return t ? t->cnt : 0;
    }

    int get_min(Tree t) {
        return t ? t->min : INF;
    }

    void update_cnt(Tree t) {
        if (t) {
            t->cnt = 1 + cnt(t->l) + cnt(t->r);
        }
    }

    void update_min(Tree t) {
        if (t) {
            t->min = min(t->value, min(get_min(t->l), get_min(t->r)));
        }
    }

    void pushup(Tree t) {
        update_cnt(t), update_min(t);
    }

    void pushdown(Tree t) {
        if (t && t->rev) {
            t->rev = false;
            swap(t->l, t->r);
            if (t->l) t->l->rev ^= 1;
            if (t->r) t->r->rev ^= 1;
        }
        if (t && t->lazy) {
            if (t->l) {
                t->l->lazy += t->lazy;
                t->l->min += t->lazy;
            }
            if (t->r) {
                t->r->lazy += t->lazy;
                t->r->min += t->lazy;
            }
            t->value += t->lazy;
            t->lazy = 0;
        }
        pushup(t);
    }
    
    void split(Tree t, int key, Tree& l, Tree& r) {
        if (!t) {
            l = r = nullptr;
            return;
        }
        pushdown(t);
        int implicit_key = cnt(t->l) + 1;
        if (key < implicit_key) {
            split(t->l, key, l, t->l), r = t;
        } else {
            split(t->r, key - implicit_key, t->r, r), l = t;
        }
        pushup(t);
    }
    
    void insert(Tree& t, int key, Tree item) {
        Tree t1, t2;
        split(t, key, t1, t2);
        merge(t1, t1, item);
        merge(t, t1, t2);
    }

    void merge(Tree& t, Tree l, Tree r) {
        pushdown(l);
        pushdown(r);
        if (!l || !r) {
            t = l ? l : r;
        } else if (l->priority > r->priority) {
            merge(l->r, l->r, r), t = l;
        } else {
            merge(r->l, l, r->l), t = r;
        }
        pushup(t);
    }
    
    void erase(Tree& t, int key) {
        Tree t1, t2, t3;
        split(t, key + 1, t1, t2);
        split(t1, key, t1, t3);
        merge(t, t1, t2);
    }

    void add(Tree t, int l, int r, T x) {
        Tree t1, t2, t3;
        split(t, l, t1, t2);
        split(t2, r - l, t2 , t3);
        t2->lazy += x;
        t2->min += x;
        merge(t2, t2, t3);
        merge(t, t1, t2);
    }

    int findmin(Tree t, int l, int r) {
        Tree t1, t2, t3;
        split(t, l, t1, t2);
        split(t2, r - l, t2, t3);
        T ret = t2->min;
        merge(t2, t2, t3);
        merge(t, t1, t2);
        return ret;
    }

    void reverse(Tree t, int l, int r) {
        if (l > r) return;
        Tree t1, t2, t3;
        split(t, l, t1, t2);
        split(t2, r - l, t2, t3);
        t2->rev ^= 1;
        merge(t2, t2, t3);
        merge(t, t1, t2);
    }

    // [l, r)の先頭がmになるように左シフトさせる。std::rotateと同じ仕様
    void rotate(Tree t, int l, int m, int r) {
        reverse(t, l, r);
        reverse(t, l, l + r - m);
        reverse(t, l + r - m, r);
    }

    void dump(Tree t) {
        if (!t) return;
        pushdown(t);
        dump(t->l);
        cout << t->value << " ";
        dump(t->r);
    }
    
public:
    ImplicitTreap(): mt(rd()), INF() {}
    explicit ImplicitTreap(T t_max)
    :ImplicitTreap(){ INF = t_max;}

    ImplicitTreap(std::vector<T> array, T t_max)
    :ImplicitTreap(t_max) {
        for(size_t idx = 0; idx < array.size(); idx++) {
            insert(idx, array[idx]);
        }
    }

    /// @brief posの前にxを追加する
    /// @param pos 
    /// @param x 
    void insert(int pos, T x) {
        insert(root, pos, new Node(x, mt(), INF));
    }

    /// @brief [l, r)の範囲にxを加算する
    /// @param l 
    /// @param r 
    /// @param x 
    void add(int l, int r, T x) {
        add(root, l, r, x);
    }

    /// @brief [l, r)の範囲の最小値を探索する
    /// @param l 
    /// @param r 
    /// @return [l, r)内の最小値
    int findmin(int l, int r) {
        return findmin(root, l, r);
    }

    /// @brief posの位置にある要素を削除する
    /// @param pos 
    void erase(int pos) {
        erase(root, pos);
    }

    /// @brief [l, r)の範囲を反転する
    /// @param l 
    /// @param r 
    void reverse(int l, int r) {
        reverse(root, l, r);
    }

    /// @brief [l, r)の範囲を、mが先頭に来るように回転する
    /// @param l 
    /// @param m 
    /// @param r 
    void rotate(int l, int m, int r) {
        rotate(root, l, m, r);
    }

    T at(int pos) {

    }

    void dump() {
        dump(root);
        cout << endl;
    }

    T operator[](int pos) {
        return findmin(pos, pos+1);
    }
};

/*封印
template<std::uint_fast64_t p=100000>
class MultiInt{
  private:
    std::vector<long long> _num;
    constexpr void shift(std::vector<long long> &a)const{
        size_t k=0;
        for(;k+1<a.size();++k){
            if(a[k]>=0){
                a[k+1] += a[k]/p;
                a[k] %= p;
            }
            else{
                long long n=(-a[k]+p-1)/p;
                a[k] += p*n;
                a[k+1] -= n;
            }
        }
        while(a[k]<=-(long long)p){
            long long n=(-a[k]+p-1)/p;
            a[k++] += p*n;
            a.push_back(-n);
        }
        while(a[k]>=(long long)p){
            a.push_back(0);
            a[k+1] += a[k]/p;
            a[k++] %= p;
        }
        while(a.back()==0) a.pop_back();
        if(a.empty()) a.push_back(0);
    }
    struct Decimal{
        MultiInt x;
        size_t k;
        Decimal cut(size_t t){if(t==x._num.size()) x=MultiInt(0);else x>>=t;k-=t;return *this;}
        Decimal(){}
        Decimal(const MultiInt &xt,const size_t &kt):x(xt),k(kt){}
        Decimal operator+(const Decimal &y)const{
            if(k<=y.k){
                MultiInt z=x;
                z<<=y.k-k;
                z+=y.x;
                return Decimal(z,y.k);
            }
            else{
                MultiInt z=y.x;
                z<<=k-y.k;
                z+=x;
                return Decimal(z,k);
            }
        }
        Decimal operator-(const Decimal &y)const{
            if(k<=y.k){
                MultiInt z=x;
                z<<=y.k-k;
                z-=y.x;
                return Decimal(z,y.k);
            }
            else{
                MultiInt z=y.x;
                z<<=k-y.k;
                z-=x;
                return Decimal(z,k);
            }
        }
        Decimal operator*(const Decimal &y)const{
            return Decimal(x*y.x,k+y.k);
        }
        Decimal operator+(const MultiInt &y)const{
            MultiInt z=y;z<<=k;z+=x;
            return Decimal(z,k);
        }
        Decimal operator-(const MultiInt &y)const{
            MultiInt z=y;z<<k;z-=x;
            return Decimal(z,k);
        }
        Decimal operator*(const MultiInt &y)const{
            return Decimal(x*y,k);
        }
        Decimal operator+=(const Decimal &y){*this = *this + y;return *this;}
        Decimal operator-=(const Decimal &y){*this = *this - y;return *this;}
        Decimal operator*=(const Decimal &y){x*=y.x;k+=y.k;return *this;}
        Decimal operator+=(const MultiInt &y){*this = *this + y;return *this;}
        Decimal operator-=(const MultiInt &y){*this = *this - y;return *this;}
        Decimal operator*=(const MultiInt &y){x*=y;return *this;}
    };
  public:
    constexpr size_t size()const{return _num.size();};
    constexpr int operator[](int k)const{return _num[k];}
    explicit constexpr MultiInt(long long num=0){
        if(num>=0){
            while(num){
                _num.push_back(num%p);
                num /= p;
            }
            if(_num.empty()) _num.push_back(0);
        }
        else{
            _num.push_back(num);
            shift(_num);
        }
    }
    explicit constexpr MultiInt(std::vector<long long> num):_num(num){shift(_num);}
    constexpr MultiInt(const MultiInt &cp):_num(cp._num){}
    constexpr MultiInt operator+(const MultiInt &x)const{
        size_t sz = std::max(_num.size(),x._num.size());
        std::vector<long long> c(sz);
        for(size_t i=0;i<sz;++i){
            if(i<_num.size()) c[i] += _num[i];
            if(i<x._num.size()) c[i] += x._num[i];
        }
        return MultiInt(c);
    }
    constexpr MultiInt operator-(const MultiInt &x)const{
        size_t sz = std::max(_num.size(),x._num.size());
        std::vector<long long> c(sz);
        for(size_t i=0;i<sz;++i){
            if(i<_num.size()) c[i] += _num[i];
            if(i<x._num.size()) c[i] -= x._num[i];
        }
        return MultiInt(c);
    }
    constexpr MultiInt operator-(void)const{
        return MultiInt(0)-*this;
    }
    constexpr MultiInt operator*(const MultiInt &x)const{
        size_t N = _num.size(),M = x._num.size();
        size_t i=0,k=0;
        if(N>=M){
            std::vector<long long> res(N+3*M),L(M);
            while(i<N){
                size_t D=i;
                for(;i-D<M;++i){
                    if(i<N) L[i-D]=_num[i];
                    else L[i-D]=0;
                }
                FFT<long long> r(L,x._num);
                for(size_t j=0;j<2*M-1;++j) res[j+k*M] += r[j];
                ++k;
            }
            return MultiInt(res);
        }
        else{
            std::vector<long long> res(M+3*N),R(N);
            while(i<M){
                size_t D=i;
                for(;i-D<N;++i){
                    if(i<M) R[i-D]=x._num[i];
                    else R[i-D]=0;
                }
                FFT<long long> r(_num,R);
                for(size_t j=0;j<2*N-1;++j) res[j+k*N] += r[j];
                ++k;
            }
            return MultiInt(res);
        }
    }
    constexpr MultiInt operator/(const MultiInt &x)const{
        MultiInt m(0),a(1),c2(2);
        size_t n=_num.size()+x._num.size();
        a<<=_num.size();
        c2<<=n;
        while(m!=a){
            m=a;
            a*=c2 - a*x;
            a>>=n;
        }
        a*=(*this);
        a>>=n;
        
        while((a+1)*x<=(*this)){
            a+=1;
        }
        return a;
    }
    constexpr MultiInt operator%(const MultiInt &x)const{return *this - (*this/x)*x;}
    constexpr MultiInt operator+(const long long t)const{
        std::vector<long long> res=_num;
        res[0]+=t;
        if(res[0]>=(long long)p) shift(res);
        return MultiInt(res);
    }
    constexpr MultiInt operator-(const long long &t)const{
        std::vector<long long> res = _num;
        res[0]-=t;
        if(res[0]<0) shift(res);
        return MultiInt(res);
    }
    constexpr MultiInt operator*(const long long &t)const{
        std::vector<long long> res=_num;
        for(size_t i=0;i<res.size();++i) res[i]*=t;
        shift(res);
        return MultiInt(res);
    }
    constexpr MultiInt operator/(const long long &t)const{return *this/MultiInt(t);}
    constexpr MultiInt operator%(const long long &t)const{return *this%MultiInt(t);}

    constexpr bool operator<(const MultiInt &x)const{
        if(_num.back()<0&&x._num.back()<0){
            if(_num.size()!=x._num.size()) return _num.size() > x._num.size();
            else for(size_t i=_num.size()-1;i>=0;--i){
                if(_num[i]!=x._num[i]) return _num[i] < x._num[i];
            }
        }
        else if(_num.back()*x._num.back()<0) return _num.back() < 0;
        else{
            if(_num.size()!=x._num.size()) return _num.size() < x._num.size();
            else for(size_t i=_num.size()-1;i>=0;--i){
                if(_num[i]!=x._num[i]) return _num[i] < x._num[i];
            }
        }
        return false;
    }
    constexpr bool operator>=(const MultiInt &x)const{return !(*this<x);}
    constexpr bool operator==(const MultiInt &x)const{
        if(_num.size()!=x._num.size()) return false;
        for(size_t i=0;i<_num.size();++i){
            if(_num[i]!=x._num[i]) return false;
        }
        return true;    
    }
    constexpr bool operator!=(const MultiInt &x)const{return _num!=x._num;}
    constexpr bool operator>(const MultiInt &x)const{return !(*this==x||*this<x);}
    constexpr MultiInt operator<<(size_t k)const{
        std::vector<long long> res(_num.size()+k);
        for(size_t i=k;i<_num.size()+k;++i) res[i]=_num[i-k];
        return MultiInt(res);
    }
    constexpr bool operator<=(const MultiInt &x)const{return !(*this>x);}
    constexpr MultiInt operator>>(size_t k)const{
        if(_num.size()<k){
            return MultiInt(0);
        }
        std::vector<long long> res(_num.size()-k);
        for(size_t i=0;i+k<_num.size();++i) res[i]=_num[i+k];
        return MultiInt(res);
    }

    constexpr MultiInt& operator+=(const MultiInt &x){
        size_t sz = std::max(_num.size(),x._num.size());
        if(_num.size()<sz) _num.resize(sz);
        for(size_t i=0;i<x._num.size();++i) _num[i]+=x._num[i];
        shift(_num);
        return *this;
    }
    constexpr MultiInt& operator-=(const MultiInt &x){
        size_t sz = std::max(_num.size(),x._num.size());
        if(_num.size()<sz) _num.resize(sz);
        for(size_t i=0;i<x._num.size();++i) _num[i]-=x._num[i];
        shift(_num);
        return *this;
    }
    constexpr MultiInt& operator*=(const MultiInt &x){*this = *this * x;return *this;}
    constexpr MultiInt& operator/=(const MultiInt &x){*this = *this / x;return *this;}
    constexpr MultiInt& operator%=(const MultiInt &x){*this = *this % x;return *this;}

    constexpr MultiInt& operator+=(const long long &a){_num[0]+=a;if(_num[0]>=(long long)p){shift(_num);}return *this;}
    constexpr MultiInt& operator-=(const long long &a){_num[0]-=a;if(_num[0]<0){shift(_num);}return *this;}
    constexpr MultiInt& operator*=(const long long &a){
        for(size_t i=0;i<_num.size();++i) _num[i]*=a;
        shift(_num);
        return *this;
    }
    constexpr MultiInt& operator/=(const long long &a){*this/=MultiInt(a);return *this;}
    constexpr MultiInt& operator<<=(size_t k){
        _num.resize(_num.size()+k);
        for(int i=_num.size()-1-k;i>=0;--i){
            _num[i+k]=_num[i];
            _num[i]=0;
        }
        return *this;
    }
    constexpr MultiInt& operator>>=(size_t k){
        if(_num.size()<k){
            _num={0};
            return *this;
        }
        for(size_t i=0;i+k<_num.size();++i) _num[i]=_num[i+k];
        while(k--) _num.pop_back();
        if(_num.empty()) _num.push_back(0);
        return *this;
    }
    explicit operator long long()const{
        long long q=1,res=0;
        for(size_t i=0;i<_num.size();++i){
            res+=_num[i]*q;q*=p;
        }
        return res;
    }
    explicit operator int()const{
        int q=1,res=0;
        for(size_t i=0;i<_num.size();++i){
            res+=_num[i]*q;q*=p;
        }
        return res;
    }
    friend constexpr MultiInt operator+(const long long a,const MultiInt &x){return x+a;}
    friend constexpr MultiInt operator-(const long long a,const MultiInt &x){return MultiInt(a)-x;}
    friend constexpr MultiInt operator*(const long long a,const MultiInt &x){return x*a;}
    friend constexpr MultiInt operator/(const long long a,const MultiInt &x){return MultiInt(a)/x;}
    friend constexpr std::ostream& operator<<(std::ostream &os,const MultiInt &x){
        size_t i=x.size();
        if(x._num.back()>=0){
            os<<x[--i];
            while(i--) printf("%05d",x[i]);
        }
        else os<<'-'<<MultiInt(-x);
        return os;
    }
    friend std::istream& operator>>(std::istream &is,MultiInt &x){
        long long val;is>>val;
        x=MultiInt(val);
        return is;
    }
};
*/