// dependancy
#include<iostream>
#include<cmath>
#include<vector>
#include<complex>
#include<stack>

/// @brief 多項式 
/// @tparam T 多項式の係数型
/// @details f(x)=a_0+a_1x+a_2x^2+...の形で表現される
/// @note fx[i]がa_iに対応する
template<typename T>
class Polynomial{
    const double PI=acos(-1);
    std::vector<T> fx;
    using Complex=std::complex<double>;

    /// @brief 高速フーリエ変換
    /// @param F 変換対象の多項式
    /// @param n Fのサイズ、2の冪乗であること
    /// @param sig 1:DFT -1:逆DFT
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
    /// @brief 逆高速フーリエ変換
    /// @param f 変換対象の多項式
    /// @param n fのサイズ、2の冪乗であること
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
    /// @brief 多項式とスカラーの積
    /// @param x スカラー
    /// @return 多項式f(x)の各係数にxをかけた多項式
    Polynomial operator*(const T &x){
        std::vector<T> hx=fx;
        for(int i=0;i<fx.size();++i){
            hx[i] = fx[i]*x;
        }
        return Polynomial(hx);
    }
    Polynomial operator+=(const Polynomial &gx){*this=*this+gx;return *this;}
    Polynomial operator-=(const Polynomial &gx){*this=*this-gx;return *this;}
    Polynomial operator*=(const Polynomial &gx){*this=(*this)*gx;return *this;}
    Polynomial operator/=(const Polynomial &gx){*this=*this/gx;return *this;}
    Polynomial operator%=(const Polynomial &gx){*this=*this%gx;return *this;}

    /// @brief 多項式の値を計算
    /// @param x 多項式に代入する値
    /// @return f(x)の値
    T operator()(const T &x)const{
        if(fx.size()==0) return 0;
        if(fx.size()==1) return fx[0]; 
        T f=fx.back();
        for(int i=fx.size()-2;i>=0;--i){
            f=f*x+fx[i];
        }
        return f;
    }
    /// @brief 多項式の積分
    /// @note 定数項は0に設定される
    void integrate(){
        fx.push_back(0);
        for(int i=fx.size()-1;i>0;--i){
            fx[i]=fx[i-1]/i;
        }
        fx[0]=0;
    }
    /// @brief 多項式の微分
    /// @note 定数項は消去される
    void differencial(){
        for(int i=0;i+1<fx.size();++i){
            fx[i]=(i+1)*fx[i+1];
        }
        fx.pop_back();
    }

    /// @brief 区間[a,b]における定積分
    /// @param a 積分区間の左端
    /// @param b 積分区間の右端
    /// @return aからbまでの定積分値
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
