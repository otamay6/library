#pragma once

#include <functional>
#include <vector>
/// @brief 区間の累積和、区間外の累積和
/// @tparam T 要素の型
/// @details 区間和をO(1)で求める
template<typename T>
class Ruiseki
{
  private:
    std::vector<T> LEFT, RIGHT;
    T d0;
    typedef std::function<T(T, T)> F;
    F f;
    F g;
    int N;

  public:
    /// @brief コンストラクタ
    /// @param a 元の配列
    /// @param f 区間和の適用関数
    /// @param g fの逆関数
    /// @param INI 零元
    Ruiseki(const std::vector<T> &a=std::vector<T>(), F f=[](T a,T b){return a+b;}, F g=[](T a,T b){return a-b;},T INI=0)
    :f(f),g(g),d0(INI)
    {
        N = a.size();
        LEFT.resize(N + 1);
        RIGHT.resize(N + 1);
        LEFT[0] = RIGHT[0] = INI;
        for(int i = 0; i < N; ++i)
        {
            LEFT[i + 1] = f(LEFT[i], a[i]);
            RIGHT[i + 1] = f(RIGHT[i], a[N - i - 1]);
        }
    }

    /// @brief [l,r)の外の累積和
    /// @param l 区間の左端(0-indexed)
    /// @param r 区間の右端(0-indexed)
    /// @return [0,l)と[r,N)の累積和
    T out(int l, int r)
    { //[l,r)の外の累積
        if(l>=r) return d0;
        if(l<0) return RIGHT[N-r];
        if(r>N) return LEFT[l];
        return F(LEFT[l], RIGHT[N - r]);
    }
    
    /// @brief [l,r)の累積和
    /// @param l 区間の左端(0-indexed)
    /// @param r 区間の右端(0-indexed)
    /// @return [l,r)の累積和
    T in(int l,int r)
    { //[l,r)内の累積
        if(l>=r) return d0;
        return g(LEFT[r],LEFT[l]);
    }
};