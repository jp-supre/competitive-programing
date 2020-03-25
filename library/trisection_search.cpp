#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
const long long INF = 1LL<<60;
const long long MOD = 1000000000+7;
#define rev(s) (string((s).rbegin(), (s).rend()))
// cout << fixed << setprecision(10) << ans << endl; 有効桁数指定
// *min_element(c + l, c + r) *max_element(c + l, c + r) 配列の中のmin-max


int main() {

    double P;
    cin >> P;

    // 目的関数(最小化したい)
    auto f = [P](double x) {
        return x + P / pow(2, x / 1.5);
    };

    // 左側
    double low = 0;
    // 右側
    double high = 1e9;

    // 500回繰り返す
    int cnt = 500;
    while (cnt--) {
        double c1 = (low * 2 + high) / 3;
        double c2 = (low + high * 2) / 3;

        // もしf(c2)のほうが良い(小さい)なら、駄目な方lowを更新する
        if (f(c1) > f(c2)) low = c1;
        else high = c2;
    }

    cout << fixed << setprecision(20) << f(low) << endl;

}