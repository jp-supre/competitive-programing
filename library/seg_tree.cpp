#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
typedef tuple<ll,ll,ll> T;
const long long INF = 1LL<<60;
const long long MOD = 1000000000+7;
#define rev(s) (string((s).rbegin(), (s).rend()))
// cout << fixed << setprecision(10) << ans << endl; 有効桁数指定
// *min_element(c + l, c + r) *max_element(c + l, c + r) 配列の中のmin-max
// int dx[8]={1,1,0,-1,-1,-1,0,1};
// int dy[8]={0,1,1,1,0,-1,-1,-1};
// int dx[4]={1,0,-1,0};
// int dy[4]={0,1,0,-1};
// ~ は、-1の時だけfalse

// INT_MAX = pow(2,31)-1;
const int MAX_N = 1 << 17;
// セグメント木を持つグローバル配列
int n, dat[2 * MAX_N - 1];
// 初期化
void init(int n_) {
// 簡単のため、要素数を2のべき乗に
    n = 1;
    while (n < n_) n *= 2;
// すべての値をINT_MAXに
    for (int i = 0; i < 2 * n - 1; i++) dat[i] = INT_MAX;
}
// k番目の値(0-indexed)をaに変更
void update(int k, int a) {
// 葉の節点
    k += n - 1;
    dat[k] = a;
// 登りながら更新
    while (k > 0) {
    k = (k - 1) / 2;
    dat[k] = min(dat[k * 2 + 1], dat[k * 2 + 2]);
    }
}
// [a, b)の最小値を求める
// 後ろのほうの引数は、計算の簡単のための引数。
// kは節点の番号、l, rはその節点が[l, r)に対応づいていることを表す。
// したがって、外からはquery(a, b, 0, 0, n)として呼ぶ。
int query(int a, int b, int k, int l, int r) {
// [a, b)と[l, r)が交差しなければ、INT_MAX
    if (r <= a || b <= l) return INT_MAX;
    // [a, b)が[l, r)を完全に含んでいれば、この節点の値
    if (a <= l && r <= b) return dat[k];
    else {
    // そうでなければ、2つの子の最小値
    int vl = query(a, b, k * 2 + 1, l, (l + r) / 2);
    int vr = query(a, b, k * 2 + 2, (l + r) / 2, r);
    return min(vl, vr);
    }
}