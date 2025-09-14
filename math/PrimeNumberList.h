#pragma once

#include<vector>
#include<utility>
#include "common.h"

class PrimeNumberList{
 private:
    std::vector<int> PN;
    std::vector<bool> isp;
 public:
    PrimeNumberList(int N){
        isp.resize(N+1,true);
        isp[0]=isp[1]=false;
        for(int i=2;i<=N;++i) if(isp[i]){
            PN.push_back(i);
            for(int j=2*i;j<=N;j+=i) isp[j]=false;
        }
    }
    int operator[](int i){return PN[i];}
    int size(){return PN.size();}
    int back(){return PN.back();}
    bool isPrime(int q){return isp[q];}
};


