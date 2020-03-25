#include <iostream>
using namespace std;

long long modinv(long long a, long long m) {
    long long b = m, u = 1, v = 0;
    while (b) {
        long long t = a / b;
        a -= t * b; swap(a, b);
        u -= t * v; swap(u, v);
    }
    u %= m; 
    if (u < 0) u += m;
    return u;
}

int main() {
    // mod. 13 での逆元を求めてみる
    for (int i = 1; i < 13; ++i) {
        cout << i << " 's inv: " << modinv(i, 13) << endl;
    }
}