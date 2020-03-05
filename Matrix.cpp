#include<iostream>
#include<vector>
using namespace std;

template<typename T> class MAT{
 private:
    int row,col;
    vector<vector<T>> _A;
    double eps = 1e-9;
    MAT set(vector<vector<T>> A){_A = A ; return *this;}
 public:
    MAT(){ }
    MAT(int n,int m=0,T x=T(0)){
        if(n<1 || m<0){cout << "err Matrix::Matrix" <<endl;exit(1);}
        row = n;
        col = m?m:n;//return E if m=0
        REP(i,row){
            vector<T> a(col,x);
            _A.push_back(a);
        }
        if(m==0) REP(i,n) _A[i][i]=1.0;
    }
    MAT(vector<vector<T>> A){row=A.size();col=A[0].size();_A=A;}
    MAT(const MAT &cp){_A=cp._A;row=cp.row;col=cp.col;}
    T* operator[] (int i){return _A[i].data();}
    MAT operator= (vector<vector<T>> x) {return set(x);}
    MAT operator+ (MAT x) const {
        if(row!=x.row || col!=x.col){
            cerr << "err Matrix::operator+" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]+x[i][j];
        return r;
    }
    MAT operator- (MAT x) const {
        if(row!=x.row || col!=x.col){
            cerr << "err Matrix::operator-" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]-x[i][j];
        return r;
    }
    MAT operator* (MAT x) const {
        if(x.col==1&&x.row==1) return x[0][0]*MAT(_A);
        if(row==1&&col==1) return _A[0][0]*x;
        if(col!=x.row){
            cerr << "err Matrix::operator*" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        MAT r(row, x.col);
        REP(i,row) REP(j,x.col) REP(k,col) r[i][j]+=_A[i][k]*x[k][j];
        return r;
    }
    MAT operator/(MAT x)const{*this = *this * x.inverse(); return *this;}
    MAT operator/ (T a)const{
        MAT r(row,col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]/a;
        return r;
    }
    MAT operator+= (MAT x) {*this = *this + x;return *this;}
    MAT operator-= (MAT x) {*this = *this - x; return *this;}
    MAT operator*=(T a){*this = a*(*this); return this;}
    MAT operator/=(MAT x){*this = *this/x;return *this;}
    MAT operator/=(T a){*this = *this/a; return *this;}
    friend MAT operator* (T n,MAT x){
        MAT r(x.row,x.col);
        REP(i,x.row) REP(j,x.col) r[i][j]=n*x[i][j];
        return r;
    }
    friend MAT operator* (MAT x,T n){
        MAT r(x.row,x.col);
        REP(i,x.row) REP(j,x.col) r[i][j]=n*x[i][j];
        return r;
    }
    explicit operator vector<vector<T>>(){return _A;}
    friend ostream &operator<<(ostream &os,const MAT &x){ REP(i,x.row) REP(j,x.col) os<<x._A[i][j]<<" \n"[j==x.col-1]; return os;}
    friend istream &operator>>(istream &is,MAT &x){REP(i,x.row) REP(j,x.col) is>>x._A[i][j];return is;}
    size_t size_row()const{return row;}
    size_t size_col()const{return col;}
    MAT transpose()const{
        MAT r(col,row);
        REP(i,col) REP(j,row) r[i][j]=_A[j][i];
        return r;
    }
    MAT inverse()const{
        T buf;
        MAT<T> inv_a(row,0);
        vector<vector<T>> a=_A;
        //row reduction
        REP(i,row){
            buf=1/a[i][i];
            REP(j,row){
                a[i][j]*=buf;
                inv_a[i][j]*=buf;
            }
            REP(j,row){
                if(i!=j){
                    buf=a[j][i];
                    REP(k,row){
                        a[j][k]-=a[i][k]*buf;
                        inv_a[j][k]-=inv_a[i][k]*buf;
                    }
                }
            }
        }
        return inv_a;
    }
    MAT Jacobi(MAT b)const{//ヤコビ法によって解を求める
        size_t sz=row;
        MAT D(sz,sz),inD(sz,sz),H(sz,sz),N(sz,sz);
        MAT c(sz,1),x(sz,1),tmp(sz,1);
        //cout<<"initialized"<<endl;
        REP(i,sz){//値の初期化、Aを対角要素とそれ以外に分解する(できてる)
            REP(j,sz){
                H[i][j] = 0;
                if(i==j){
                    D[i][j] = _A[i][j];
                    inD[i][j] = 1/_A[i][j];
                    N[i][j]=0;
                }
                else if(i!=j){
                    D[i][j] = 0;
                    inD[i][j] = 0;
                    N[i][j]=_A[i][j];
                }
            }
            c[i][0] = 0;
            x[i][0] = 1;
        }
        c=inD*b;
        H=inD*N;
        while(1){//反復法ステップ1→2→1...
            tmp=x;
            x=c-H*x;
            T r=T(0);
            for(int i=0;i<row;++i){
                r+=(x[i][0]-tmp[i][0])*(x[i][0]-tmp[i][0]);
            }
            if(r<eps) break;
        }
        return x;
    }
    MAT Gauss(MAT b)const{//ガウス・ザイデル法によって解を求める
        MAT<T> DL(row),U(row),inDL(row),H(row),c(row,1),x(row,1),tmp(row,1);
        for(int i=0;i<row;i++){
            for(int j=0;j<col;j++){
                H[i][j] = 0;
                if(i>=j){
                    DL[i][j] = _A[i][j];
                    U[i][j] = 0;
                }
                else{
                    DL[i][j] = 0;
                    U[i][j] = _A[i][j];
                }
            }
            x[i][0] = 1;
        }
        inDL=DL.inverse();
        c=inDL*b;
        H=inDL*U;
        int n=0;
        while(1){
            tmp=x;
            x=c-H*x;
            T r = T(0);
            for(int i=0;i<row;++i){
                r+=(x[i][0]-tmp[i][0])*(x[i][0]-tmp[i][0]);
            }
            n++;
            if(r<eps) break;
        }
        return x;
    }
    int rank()const{// O( n^3 )
        vector<vector<T>> A=_A;
        const int n = row, m = col;
        int r = 0;
        for(int i = 0; r < n && i < m; ++i) {
            int pivot = r;
            for(int j = r+1; j < n; ++j) if(fabs(A[j][i]) > fabs(A[pivot][i])) pivot = j;
            swap(A[pivot], A[r]);
            if(fabs(A[r][i]) < eps) continue;
            for (int k = m-1; k >= i; --k) A[r][k] /= A[r][i];
            rep(j,r+1,n) rep(k,i,m) A[j][k] -= A[r][k] * A[j][i];
            ++r;
        }
        return r;
    }
};

template<typename T>
struct CSR{
    std::map<std::pair<int,int>,T> List;
    std::vector<T> A,IA,JA,nIA;
    int H,W;
    CSR(int H,int W):H(H),W(W){}
    CSR(const CSR &X){*this=X;}
    void add_val(int row,int col,T val){
        List[std::make_pair(row,col)] += val;
    }
    CSR Unit()const{
        CSR res(H,W);
        for(int i=0;i<H;++i){
            res.add_val(i,i,1);
        }
        res.compress();
        return res;
    }
    void compress(){
        int N=List.size();
        A.resize(N);
        IA.resize(N);
        JA.resize(N);
        int ind=0;
        for(auto L:List){
            std::pair<int,int> pos=L.first;
            A[ind]=L.second;
            IA[ind]=pos.first;
            JA[ind]=pos.second;
            ind++;
        }
        nIA.resize(H+1,List.size());
        for(int i=0,r=0;i<IA.size();++i){
            while(r<=IA[i]){
                nIA[r]=i;
                r++;
            }
        }
    }
    CSR operator*(const CSR &X)const{
        if(W!=X.H){
            std::cerr << "err Matrix::operator*" <<std::endl;
            std::cerr << "  not equal matrix size" <<std::endl;
            exit(0);
        }
        CSR res(H,X.W);
        for(int i=0;i<H;++i){
            for(int j=nIA[i];j<nIA[i+1];++j){
                int k=JA[j];
                for(int t=X.nIA[k];t<X.nIA[k+1];++t){
                    res.add_val(i,X.JA[t],A[j]*X.A[t]);
                }
            }
        }
        res.compress();
        return res;
    }
    friend std::ostream& operator<<(std::ostream &os,CSR &X){
        for(int i=0;i<X.H;++i) for(int j=0;j<X.W;++j){
            std::pair<int,int> p=std::make_pair(i,j);
            if(X.List.count(p)) os<<X.List[p];
            else os<<0;
            os<<" \n"[j+1==X.W&&i+1<X.H];
        }
        return os;
    }
};