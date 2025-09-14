#pragma once

#include <iostream>
#include <string>
#include <cstring>

namespace strlib {
    /// @brief 文字列をa桁で分割
    /// @param n 分割対象の文字列
    /// @param a 分割する桁数
    /// @return (左側の文字列,右側の文字列)、aが不正な場合(-1,-1)
    inline std::pair<long long, long long> splitint(const std::string &n,int a){
        int Len=n.length();
        if(a<0||Len<a) return {-1,-1};
        std::string left,right;
        if(a!=0) left=n.substr(0,a);
        if(a!=Len) right=n.substr(a);
        return std::make_pair(stoll(left),stoll(right));
    }

    /// @brief 回文判定
    /// @param s 判定対象の文字列
    /// @return 回文ならtrue
    inline bool is_kaibun(const std::string &s){//O(|S|)
        size_t N = s.size();
        for(size_t i = 0; i < N/2; ++i){
            if(s[i] != s[N-i-1]) return false;
        }
        return true;
    }

    inline void printYesNo(bool yes){
        std::cout << (yes ? "Yes\n" : "No\n");
    }
};