#include <cassert>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <climits>

enum class OverFlowBehave{
    ROLL,
    ERROR,  // SafeInt自体は自発的にエラー通知をしない。SATURATEでの代用を推奨
    SATURATE
};

template<int MIN, int MAX, OverFlowBehave overflow_behave = OverFlowBehave::ROLL>
class SafeInt{
    static_assert(INT_MIN < MIN && MAX < INT_MAX && INT_MAX - MAX - 1 <= -MIN, "SafeIntは本来のint型より明確に小さい範囲で使用してください");
private:
    int value;
public:
    SafeInt(int v = 0) {
        if(OverFlowBehave::SATURATE == overflow_behave){
            if(v < MIN) v = MIN;
            if(v > MAX) v = MAX;
        } else if(OverFlowBehave::ROLL == overflow_behave){
            constexpr int cycle = MAX - MIN + 1;
            if(v < MIN){
                // diffの算出自体がオーバーフローしないか
                // MAXと逆のやり方でいい
                if(v < 0 && INT_MAX + v < MIN){
                    v += cycle * (-v / cycle);
                    if(INT_MAX + v < MIN)
                        v += cycle;
                }
                int diff = MIN - v;
                diff %= cycle;
                v = MAX - diff + 1;
            }
            if(v > MAX){
                // diffの算出自体がオーバーフローしないか
                if(v > 0 && INT_MAX - v < - MAX){
                    // MAXが負は自明なのでvが正ならばMAXより大きい
                    v %= cycle;
                    // 尚オーバーフローしそうなら更に引く
                    if(INT_MAX - v < - MAX)
                        v -= cycle;
                }
                int diff = v - MAX;
                diff %= cycle;
                v = MIN + diff - 1;
            }
        } else if(OverFlowBehave::ERROR == overflow_behave){
            v = INT_MAX;
        }
        value = v;
    }
    /// @brief オーバーフローしているか判定
    /// @return オーバーフロー状態かどうか
    bool overflowed() const {
        if constexpr (overflow_behave == OverFlowBehave::ERROR){
            return value == INT_MAX;
        }
        return false;
    }
    // 四則演算
    friend SafeInt operator +(const SafeInt &a, const SafeInt &b) {
        if(a.overflowed() || b.overflowed()) return INT_MAX;
        int result = a.value + b.value;
        if(a.value >= 0 && b.value >= 0) {
            if(MAX - a.value < b.value){
                if constexpr (overflow_behave == OverFlowBehave::ROLL) {
                    result = MIN + (b.value - (MAX - a.value) - 1);
                } else if(overflow_behave == OverFlowBehave::SATURATE){
                    result = MAX;
                } else if(overflow_behave == OverFlowBehave::ERROR){
                    result = INT_MAX;
                }
            }
        }
        else if(a.value < 0 && b.value < 0) {
            if(MIN - a.value > b.value){
                if constexpr (overflow_behave == OverFlowBehave::ROLL) {
                    result = MAX + (b.value - (MIN - a.value) + 1);
                } else if(overflow_behave == OverFlowBehave::SATURATE){
                    result = MIN;
                } else if(overflow_behave == OverFlowBehave::ERROR){
                    result = INT_MAX;
                }
            }
        }
        return SafeInt(result);
    }
    SafeInt operator -() const {
        return SafeInt(0) - (*this);
    }
    friend SafeInt operator -(const SafeInt &a, const SafeInt &b) {
        if(a.overflowed() || b.overflowed()) return INT_MAX;
        int result = a.value - b.value;
        if(a.value >= 0 && b.value < 0) {
            if(MAX + b.value < a.value){
                if constexpr (overflow_behave == OverFlowBehave::ROLL) {
                    result = MIN + (-b.value - (MAX - a.value) - 1);
                } else if(overflow_behave == OverFlowBehave::SATURATE){
                    result = MAX;
                } else if(overflow_behave == OverFlowBehave::ERROR){
                    result = INT_MAX;
                }
            }
        }
        else if(a.value < 0 && b.value >= 0) {
            if(MIN + b.value > a.value){
                if constexpr (overflow_behave == OverFlowBehave::ROLL) {
                    result = MAX - (b.value - (MIN - a.value) + 1);
                } else if(overflow_behave == OverFlowBehave::SATURATE){
                    result = MIN;
                } else if(overflow_behave == OverFlowBehave::ERROR){
                    result = INT_MAX;
                }
            }
        }
        return SafeInt(result);
    }
    friend SafeInt operator *(const SafeInt &a, const SafeInt &b) {
        if(a.overflowed() || b.overflowed()) return INT_MAX;
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
                if(result.overflowed()){
                    break;
                }
                if(!first && overflow_behave == OverFlowBehave::SATURATE){
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
        if(a.overflowed() || b.overflowed()) return INT_MAX;
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
    friend SafeInt operator %(const SafeInt &a, const SafeInt &b){
        if(a.overflowed() || b.overflowed()) return INT_MAX;
        // 剰余はbが正であることを仮定していい
        assert(b.value > 0);
        // 最小値が0以下でないと剰余を取る意味がない
        assert(MIN <= 0);
        return SafeInt(a.value % b.value);
    }
    // 比較演算
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

    SafeInt& operator++(){
        if(this->overflowed()) return *this;
        *this += 1;
        return *this;
    }
    SafeInt operator ++(int){
        if(this->overflowed()) return *this;
        SafeInt temp = *this;
        ++*this;
        return temp;
    }

    SafeInt& operator--(){
        if(this->overflowed()) return *this;
        *this -= 1;
        return *this;
    }
    SafeInt operator --(int){
        if(this->overflowed()) return *this;
        SafeInt temp = *this;
        --*this;
        return temp;
    }
    SafeInt operator +=(const SafeInt &x){
        *this = *this + x;
        return *this;
    }
    SafeInt operator -=(const SafeInt &x){
        *this = *this - x;
        return *this;
    }
    SafeInt operator *=(const SafeInt &x){
        *this = *this * x;
        return *this;
    }
    SafeInt operator /=(const SafeInt &x){
        *this = *this / x;
        return *this;
    }
    SafeInt operator %=(const SafeInt &x){
        *this = *this % x;
        return *this;
    }
    // 整数変換
    explicit operator int(){
        return this->value;
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
    // 出力
    friend std::ostream& operator<<(std::ostream &os, const SafeInt &a) {
        os << a.value;
        return os;
    }
    friend std::istream& operator>>(std::istream &is, SafeInt &a) {
        int v;
        is >> v;
        assert(MIN <= v && v <= MAX);
        a = SafeInt(v);
        return is;
    }
};
