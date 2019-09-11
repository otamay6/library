#include<iostream>
typedef long long ll;
constexpr ll mod=1e9+7;
class mint{
  private:
    ll _num, _mod;
    mint set(ll num)
    {
        _num = num;
        if (_num >= 0)
            _num %= _mod;
        else
            _num += (1 - (_num + 1) / _mod) * _mod;
        return *this;
    }
    ll _mpow(ll x, ll n)
    { 
        ll ans = 1;
        while (n != 0)
        {
            if (n & 1)
                ans = ans * x % _mod;
            x = x * x % _mod;
            n = n >> 1;
        }
        return ans;
    }
    ll imod(ll n) { return _mpow(n, _mod - 2); }

  public:
    mint()
    {
        _num = 0;
        _mod = mod;
    }
    mint(ll num)
    {
        _mod = mod;
        _num = (num + (1LL << 25) * mod) % mod;
    }
    mint(ll num, ll M)
    {
        _mod = M;
        _num = (num + (1LL << 25) * mod) % _mod;
    }
    mint(const mint &cp)
    {
        _num = cp._num;
        _mod = cp._mod;
    }
    mint operator=(const ll x) { return set(x); }
    mint operator+(const ll x) { return mint(_num + (x % _mod), _mod); }
    mint operator-(const ll x) { return mint(_num - (x % _mod), _mod); }
    mint operator*(const ll x) { return mint(_num * (x % _mod), _mod); }
    mint operator/(ll x) { return mint(_num * imod(x), _mod); }
    mint operator+=(const ll x) { return set(_num + (x % _mod)); }
    mint operator-=(const ll x) { return set(_num - (x % _mod)); }
    mint operator*=(const ll x) { return set(_num * (x % _mod)); }
    mint operator/=(ll x) { return set(_num * imod(x)); }
    mint operator+(const mint &x) { return mint(_num + x._num, _mod); }
    mint operator-(const mint &x) { return mint(_num - x._num, _mod); }
    mint operator*(const mint &x) { return mint(_num * x._num, _mod); }
    mint operator/(mint x) { return mint(_num * imod(x._num), _mod); }
    mint operator+=(const mint &x) { return set(_num + x._num); }
    mint operator-=(const mint &x) { return set(_num - x._num); }
    mint operator*=(const mint &x) { return set(_num * x._num); }
    mint operator/=(mint x) { return set(_num * imod(x._num)); }

    bool operator<(const mint &x) const { return _num < x._num; }
    bool operator==(const mint &x) const { return _num == x._num; }
    bool operator>(const mint &x) const { return _num > x._num; }

    friend mint operator+(ll x, const mint &m) { return mint(m._num + (x % m._mod), m._mod); }
    friend mint operator-(ll x, const mint &m) { return mint((x % m._mod) - m._num, m._mod); }
    friend mint operator*(ll x, const mint &m) { return mint(m._num * (x % m._mod), m._mod); }
    friend mint operator/(ll x, mint m) { return mint(m.imod(m._num) * x, m._mod); }

    explicit operator ll() { return _num; }
    explicit operator int() { return (int)_num; }

    friend std::ostream &operator<<(std::ostream &os, const mint &x)
    {
        os << x._num;
        return os;
    }
    friend std::istream &operator>>(std::istream &is, mint &x)
    {
        ll val;
        is >> val;
        x.set(val);
        return is;
    }
};
struct rational
{
    long long p, q;
    long long gcd(ll x,ll y){return y?gcd(y,x%y):x;}
    void normalize()
    { // keep q positive
        if (q < 0)
            p *= -1, q *= -1;
        long long d = gcd(p < 0 ? -p : p, q);
        if (d == 0)
            p = 0, q = 1;
        else
            p /= d, q /= d;
    }
    rational(long long p, long long q = 1) : p(p), q(q)
    {
        normalize();
    }
    rational &operator+=(const rational &a)
    {
        p = a.q * p + a.p * q;
        q = a.q * q;
        normalize();
        return *this;
    }
    rational &operator-=(const rational &a)
    {
        p = a.q * p - a.p * q;
        q = a.q * q;
        normalize();
        return *this;
    }
    rational &operator*=(const rational &a)
    {
        p *= a.p;
        q *= a.q;
        normalize();
        return *this;
    }
    rational &operator/=(const rational &a)
    {
        p *= a.q;
        q *= a.p;
        normalize();
        return *this;
    }
    rational &operator-()
    {
        p *= -1;
        return *this;
    }
    friend rational operator+(const rational &a, const rational &b) { return rational(a) += b; }
    friend rational operator*(const rational &a, const rational &b) { return rational(a) *= b; }
    friend rational operator-(const rational &a, const rational &b) { return rational(a) -= b; }
    friend rational operator/(const rational &a, const rational &b) { return rational(a) /= b; }
    friend bool operator<(const rational &a, const rational &b)
    { // avoid overflow
        return (long double)a.p * b.q < (long double)a.q * b.p;
    }
    friend bool operator<=(const rational &a, const rational &b) { return !(b < a); }
    friend bool operator>(const rational &a, const rational &b) { return b < a; }
    friend bool operator>=(const rational &a, const rational &b) { return !(a < b); }
    friend bool operator==(const rational &a, const rational &b) { return !(a < b) && !(b < a); }
    friend bool operator!=(const rational &a, const rational &b) { return (a < b) || (b < a); }
    friend std::ostream &operator<<(std::ostream &os, const rational &x)
    {
        printf("%.16f", (double)x.p / (double)x.q);
        return os;
    }
    friend std::istream &operator>>(std::istream &is, rational &x)
    {
        is >> x.p >> x.q;
        x.normalize();
        return is;
    }
};