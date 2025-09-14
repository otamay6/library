#pragma once

#include "common.h"
#include "PrimeNumberList.h"

class DivisorList
{ //素因数分解をしてから約数列挙、分解結果はＰ(底,指数)でpfacにまとめている
private:
    std::vector<mathlib::ll> F;
    std::vector<std::pair<int,int>> pfactorize;

public:
    DivisorList(mathlib::ll N){
        for(mathlib::ll i=1; i*i<=N; i++){
            if(N%i == 0){
                F.push_back(i);
                if (i*i != N) F.push_back(N/i);
            }
        }
        std::sort(begin(F), end(F));
        PrimeNumberList p(mathlib::sqrt_floor(N) + 1);
        for(int i=0;i< p.size();++i){
            pfactorize.push_back(std::make_pair(p[i], 0));
            while(N%p[i] == 0){
                N /= p[i];
                pfactorize.back().second++;
            }
            if (pfactorize.back().second == 0)
                pfactorize.pop_back();
        }
        if (N > 1)
            pfactorize.push_back(std::make_pair(N, 1));
    }
    int size() { return F.size(); }
    const std::vector<std::pair<int,int>> &pfac() { return pfactorize; }
    const mathlib::ll &operator[](int k) const { return F[k]; }
};
