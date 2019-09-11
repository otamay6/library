#include<iostream>
#include<vector>
#include<math.h>
#include<utility>
#define REP(i,n) for(int i=0;i<n;++i)
using namespace std;
template <typename T>
class MAT
{
  private:
    int row, col;
    vector<vector<T>> _A;
    MAT set(vector<vector<T>> A)
    {
        _A = A;
        return *this;
    }

  public:
    MAT() {}
    MAT(int n, int m)
    {
        if (n < 1 || m < 0)
        {
            cout << "err Matrix::Matrix" << endl;
            exit(1);
        }
        row = n;
        col = m ? m : n; //m=0のとき単位行列を作る
        REP(i, row)
        {
            vector<T> a(col, 0.0);
            _A.push_back(a);
        }
        //  値の初期化
        if (m == 0)
            REP(i, n)
            _A[i][i] = 1.0;
    }
    MAT(const MAT &cp)
    {
        _A = cp._A;
        row = cp.row;
        col = cp.col;
    }
    T *operator[](int i) { return _A[i].data(); }
    MAT operator=(vector<vector<T>> x) { return set(x); }
    MAT operator+(MAT x)
    {
        if (row != x.row || col != x.col)
        {
            cout << "err Matrix::operator+" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] + x[i][j];
        return r;
    }
    MAT operator-(MAT x)
    {
        if (row != x.row || col != x.col)
        {
            cout << "err Matrix::operator-" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] - x[i][j];
        return r;
    }
    MAT operator*(MAT x)
    {
        if (col != x.row)
        {
            cout << "err Matrix::operator*" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, x.col);
        REP(i, row)
        REP(j, x.col) REP(k, col) r[i][j] += _A[i][k] * x[k][j];
        return r;
    }
    MAT operator/(T a)
    {
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] / a;
        return r;
    }
    MAT operator^(ll n)
    {
        if (row != col)
        {
            cout << "err Matrix::operator^" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, 0), A = *this;
        while (n)
        {
            if (n & 1)
                r *= A;
            A *= A;
            n >>= 1;
        }
        return r;
    }
    MAT operator+=(MAT x)
    {
        if (row != x.row || col != x.col)
        {
            cout << "err Matrix::operator+=" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] + x[i][j];
        return set(r._A);
    }
    MAT operator-=(MAT x)
    {
        if (row != x.row || col != x.col)
        {
            cout << "err Matrix::operator-=" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] - x[i][j];
        return set(r._A);
    }
    MAT operator*=(MAT x)
    {
        if (col != x.row)
        {
            cout << "err Matrix::operator*" << endl;
            cout << "  not equal matrix size" << endl;
            exit(0);
        }
        MAT r(row, x.col);
        REP(i, row)
        REP(j, x.col) REP(k, col) r[i][j] += _A[i][k] * x[k][j];
        return set(r._A);
    }
    MAT operator/=(T a)
    {
        MAT r(row, col);
        REP(i, row)
        REP(j, col) r[i][j] = _A[i][j] / a;
        return r;
    }

    friend MAT operator*(T n, MAT x)
    {
        MAT r(x.row, x.col);
        REP(i, x.row)
        REP(j, x.col) r[i][j] = n * x[i][j];
        return r;
    }
    friend MAT operator*(MAT x, T n)
    {
        MAT r(x.row, x.col);
        REP(i, x.row)
        REP(j, x.col) r[i][j] = n * x[i][j];
        return r;
    }
    explicit operator vector<vector<T>>() { return _A; }
    friend ostream &operator<<(ostream &os, const MAT &x)
    {
        REP(i, x.row)
            REP(j, x.col) os << x._A[i][j] << " \n"[j == x.col - 1];
        return os;
    }
    friend istream &operator>>(istream &is, MAT &x)
    {
        REP(i, x.row)
            REP(j, x.col) is >> x._A[i][j];
        return is;
    }
    int size_row() { return row; }
    int size_col() { return col; }
    MAT transpose()
    {
        MAT r(col, row);
        REP(i, col)
        REP(j, row) r[i][j] = _A[j][i];
        return r;
    }
    MAT inverse()
    {
        T buf;
        MAT<T> inv_a(row, 0);
        vector<vector<T>> a = _A;
        //掃き出し法
        REP(i, row)
        {
            buf = 1 / a[i][i];
            REP(j, row)
            {
                a[i][j] *= buf;
                inv_a[i][j] *= buf;
            }
            REP(j, row)
            {
                if (i != j)
                {
                    buf = a[j][i];
                    REP(k, row)
                    {
                        a[j][k] -= a[i][k] * buf;
                        inv_a[j][k] -= inv_a[i][k] * buf;
                    }
                }
            }
        }
        return inv_a;
    }
    // O( n^3 ).
    int rank()
    {
        vector<vector<T>> A = _A;
        const int n = row, m = col;
        int r = 0;
        for (int i = 0; r < n && i < m; ++i)
        {
            int pivot = r;
            for (int j = r + 1; j < n; ++j)
                if (fabs(A[j][i]) > fabs(A[pivot][i]))
                    pivot = j;
            swap(A[pivot], A[r]);
            if (fabs(A[r][i]) < eps)
                continue;
            for (int k = m - 1; k >= i; --k)
                A[r][k] /= A[r][i];
            rep(j, r + 1, n) rep(k, i, m) A[j][k] -= A[r][k] * A[j][i];
            ++r;
        }
        return r;
    }
};