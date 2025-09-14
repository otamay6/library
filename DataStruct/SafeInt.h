#include <iostream>
#include <stdexcept>
#include <vector>
#include <climits>

template<int MIN, int MAX, bool rolling = false>
class SafeInt{
    static_assert(INT_MIN < MIN && MAX < INT_MAX, "SafeIntは本来のint型より小さい範囲で使用してください");
private:
    int value;
public:
    explicit SafeInt(int v = 0) {
        if(v < MIN) v = MIN;
        if(v > MAX) v = MAX;
        value = v;
    }
    friend SafeInt operator +(const SafeInt &a, const SafeInt &b) {
        int result = a.value + b.value;
        if(a.value >= 0 && b.value >= 0) {
            if(MAX - a.value < b.value){
                if constexpr (rolling) {
                    result = MIN + (b.value - (MAX - a.value) - 1);
                } else {
                    result = MAX;
                }
            }
        }
        else if(a.value < 0 && b.value < 0) {
            if(MIN - a.value > b.value){
                if constexpr (rolling) {
                    result = MAX + (b.value - (MIN - a.value) + 1);
                } else {
                    result = MIN;
                }
            }
        }
        return SafeInt(result);
    }
    SafeInt operator -() const {
        return SafeInt(0) - (*this);
    }
    friend SafeInt operator -(const SafeInt &a, const SafeInt &b) {
        int result = a.value - b.value;
        if(a.value >= 0 && b.value < 0) {
            if(MAX + b.value < a.value){
                if constexpr (rolling) {
                    result = MIN + (-b.value - (MAX - a.value) - 1);
                } else {
                    result = MAX;
                }
            }
        }
        else if(a.value < 0 && b.value >= 0) {
            if(MIN + b.value > a.value){
                if constexpr (rolling) {
                    result = MAX - (b.value - (MIN - a.value) + 1);
                } else {
                    result = MIN;
                }
            }
        }
        return SafeInt(result);
    }
    friend SafeInt operator *(const SafeInt &a, const SafeInt &b) {
        // 二分累乗法でaをb回足する
        if(abs(a.value) < abs(b.value)) return b * a;
        if(b.value == 0) return SafeInt(0);
        int n = abs(b.value);
        SafeInt result, x = b.value >= 0 ? a : -a;
        bool first = true;
        while(n) {
            if(n & 1) {
                if(first){
                    result = x;
                    first = false;
                }
                else {
                    result = result + x;
                }
                if(!first && !rolling){
                    if(result.value == MAX && x.value >= 0) break;
                    if(result.value == MIN && x.value <= 0) break;
                }
            }
            n >>= 1;
            x = x + x;
        }
        return result;
    }
    friend SafeInt operator /(const SafeInt &a, const SafeInt &b) {
        if(b.value == 0){
            if(a.value > 0) return SafeInt(MAX);
            if(a.value < 0) return SafeInt(MIN);
            return SafeInt(0);
        }
        // 範囲が逸脱した場合roolingせずにMAXかMINを返す
        if(a.value > 0 && b.value < 0) {
            if(a.value / b.value < MIN){
                return SafeInt(MIN);
            }
        }
        if(a.value < 0 && b.value < 0){
            if(a.value / b.value > MAX){
                return SafeInt(MAX);
            }
        }
        return SafeInt(a.value / b.value);
    }
    friend bool operator ==(const SafeInt &a, const SafeInt &b) {
        return a.value == b.value;
    }
    friend bool operator !=(const SafeInt &a, const SafeInt &b) {
        return a.value != b.value;
    }
    friend bool operator <(const SafeInt &a, const SafeInt &b) {
        return a.value < b.value;
    }
    friend bool operator >(const SafeInt &a, const SafeInt &b) {
        return a.value > b.value;
    }
    friend bool operator <=(const SafeInt &a, const SafeInt &b) {
        return a.value <= b.value;
    }
    friend bool operator >=(const SafeInt &a, const SafeInt &b) {
        return a.value >= b.value;
    }
    friend std::ostream& operator<<(std::ostream &os, const SafeInt &a) {
        os << a.value;
        return os;
    }
    bool operator ==(const int &b) const {
        return value == b;
    }
    bool operator !=(const int &b) const {
        return value != b;
    }
    bool operator <(const int &b) const {
        return value < b;
    }
    bool operator >(const int &b) const {
        return value > b;
    }
    bool operator <=(const int &b) const {
        return value <= b;
    }
    bool operator >=(const int &b) const {
        return value >= b;
    }
    friend std::istream& operator>>(std::istream &is, SafeInt &a) {
        int v;
        is >> v;
        a = SafeInt(v);
        return is;
    }
};
