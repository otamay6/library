// dependancies
#include<iostream>
#include<vector>
#include<map>
using namespace std;
#define REP(i,n) for(int i=0;i<(n);++i)
#define rep(i,a,b) for(int i=int(a);i<int(b);++i)

/// @brief 行列
/// @tparam T 行列の要素型
/// @details 行列の各要素はT型で表現される
template<typename T> class Matrix{
 private:
    int row,col;
    vector<vector<T>> _A;
    double eps = 1e-9;
    Matrix set(vector<vector<T>> A){_A = A ; return *this;}
 public:
    Matrix(){ }
    /// @brief n行m列の行列を生成、m=0なら単位行列
    /// @param n 行数
    /// @param m 列数
    /// @param x 行列の初期値
    explicit Matrix(int n,int m=0,T x=T(0)){
        if(n<1 || m<0){cout << "err Matrix::Matrix" <<endl;exit(1);}
        row = n;
        col = m?m:n;//return E if m=0
        REP(i,row){
            vector<T> a(col,x);
            _A.push_back(a);
        }
        if(m==0) REP(i,n) _A[i][i]=1.0;
    }
    /// @brief 行列を生成
    /// @param A 行列の初期値
    explicit Matrix(const vector<vector<T>> &A){row=A.size();col=A[0].size();_A=A;}
    Matrix(const Matrix &cp){_A=cp._A;row=cp.row;col=cp.col;}
    /// @brief 行列の各要素にアクセス
    /// @param i 行インデックス
    /// @return i行目のvector
    T* operator[] (int i){return _A[i].data();}
    Matrix operator= (vector<vector<T>> x) {return set(x);}
    Matrix operator+ (const Matrix &x) const {
        if(row!=x.row || col!=x.col){
            cerr << "err Matrix::operator+" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        Matrix r(row, col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]+x[i][j];
        return r;
    }
    Matrix operator- (const Matrix &x) const {
        if(row!=x.row || col!=x.col){
            cerr << "err Matrix::operator-" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        Matrix r(row, col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]-x[i][j];
        return r;
    }
    Matrix operator* (const Matrix &x) const {
        if(x.col==1&&x.row==1) return x[0][0]*Matrix(_A);
        if(row==1&&col==1) return _A[0][0]*x;
        if(col!=x.row){
            cerr << "err Matrix::operator*" <<endl;
            cerr << "  not equal matrix size" <<endl;
            exit(0);
        }
        Matrix r(row, x.col);
        REP(i,row) REP(j,x.col) REP(k,col) r[i][j]+=_A[i][k]*x[k][j];
        return r;
    }
    Matrix operator/(const Matrix &x)const{*this = *this * x.inverse(); return *this;}
    Matrix operator/ (const T &a)const{
        Matrix r(row,col);
        REP(i,row) REP(j,col) r[i][j]=_A[i][j]/a;
        return r;
    }
    Matrix operator+= (const Matrix &x) {*this = *this + x;return *this;}
    Matrix operator-= (const Matrix &x) {*this = *this - x; return *this;}
    Matrix operator*=(const T &a){*this = a*(*this); return this;}
    Matrix operator/=(const Matrix &x){*this = *this/x;return *this;}
    Matrix operator/=(const T &a){*this = *this/a; return *this;}
    friend Matrix operator* (const T &n,const Matrix &x){
        Matrix r(x.row,x.col);
        REP(i,x.row) REP(j,x.col) r[i][j]=n*x[i][j];
        return r;
    }
    friend Matrix operator* (const Matrix &x,const T &n){
        Matrix r(x.row,x.col);
        REP(i,x.row) REP(j,x.col) r[i][j]=n*x[i][j];
        return r;
    }
    explicit operator std::vector<vector<T>>(){return _A;}
    friend ostream &operator<<(ostream &os,const Matrix &x){ REP(i,x.row) REP(j,x.col) os<<x._A[i][j]<<" \n"[j==x.col-1]; return os;}
    friend istream &operator>>(istream &is,Matrix &x){REP(i,x.row) REP(j,x.col) is>>x._A[i][j];return is;}
    size_t size_row()const{return row;}
    size_t size_col()const{return col;}
    /// @brief 転置行列を返す
    Matrix transpose()const{
        Matrix r(col,row);
        REP(i,col) REP(j,row) r[i][j]=_A[j][i];
        return r;
    }
    /// @brief 逆行列を返す
    Matrix inverse()const{
        T buf;
        Matrix<T> inv_a(row,0);
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
    /// @brief ヤコビ法によって連立一次方程式Ax=bを解く
    /// @param b 連立方程式の定数ベクトル
    /// @return 連立方程式の解ベクトル
    /// @note 計算量はO(n^2)程度
    Matrix Jacobi(const Matrix &b)const{//ヤコビ法によって解を求める
        size_t sz=row;
        Matrix D(sz,sz),inD(sz,sz),H(sz,sz),N(sz,sz);
        Matrix c(sz,1),x(sz,1),tmp(sz,1);
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

    /// @brief 連立一次方程式をガウス・ザイデル法で解く
    /// @param b 連立方程式の定数ベクトル
    /// @return 連立方程式の解ベクトル
    /// @details 収束しない場合がある
    /// @note 計算量はO(n^2)程度
    Matrix Gauss(const Matrix &b)const{//ガウス・ザイデル法によって解を求める
        Matrix<T> DL(row),U(row),inDL(row),H(row),c(row,1),x(row,1),tmp(row,1);
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

    /// @brief 行列の階数を求める
    /// @return 行列の階数
    /// @note 計算量はO(n^3)程度
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

/// @brief 疎行列
/// @tparam T 行列の要素型
/// @details 行列の各要素はT型で表現される
template<typename T>
struct CSRMatrix{
    std::map<std::pair<int,int>,T> List;
    std::vector<T> A,IA,JA,nIA;
    int H,W;
    CSRMatrix(int H,int W):H(H),W(W){}
    CSRMatrix(const CSRMatrix &X){*this=X;}
    void add_val(int row,int col,T val){
        List[std::make_pair(row,col)] += val;
    }
    /// @brief 単位行列を返す
    CSRMatrix Unit()const{
        CSRMatrix res(H,W);
        for(int i=0;i<H;++i){
            res.add_val(i,i,1);
        }
        res.compress();
        return res;
    }
    /// @brief CSR形式に変換する
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
    /// @brief 行列積を返す
    /// @param X かける行列
    /// @return 積
    CSRMatrix operator*(const CSRMatrix &X)const{
        if(W!=X.H){
            std::cerr << "err Matrix::operator*" <<std::endl;
            std::cerr << "  not equal matrix size" <<std::endl;
            exit(0);
        }
        CSRMatrix res(H,X.W);
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
    friend std::ostream& operator<<(std::ostream &os,CSRMatrix &X){
        for(int i=0;i<X.H;++i) for(int j=0;j<X.W;++j){
            std::pair<int,int> p=std::make_pair(i,j);
            if(X.List.count(p)) os<<X.List[p];
            else os<<0;
            os<<" \n"[j+1==X.W&&i+1<X.H];
        }
        return os;
    }
};