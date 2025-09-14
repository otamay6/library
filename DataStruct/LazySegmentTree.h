#pragma once
#include <vector>
#include <functional>
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
