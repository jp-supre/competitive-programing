#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
typedef pair<int,P> Pp;
const long long INF = 1LL<<60;
const long long MOD = 1000000000+7;
#define rev(s) (string((s).rbegin(), (s).rend()))
// cout << fixed << setprecision(10) << ans << endl; 有効桁数指定
// *min_element(c + l, c + r) *max_element(c + l, c + r) 配列の中のmin-max
 
// UnionFind tree(n);  と宣言
 
struct UnionFind {
    vector<int> par; // par[i]:iの親の番号　(例) par[3] = 2 : 3の親が2
 
    UnionFind(int N) : par(N) { //最初は全てが根であるとして初期化
        for(int i = 0; i < N; i++) par[i] = i;
    }
 
    int root(int x) { // データxが属する木の根を再帰で得る：root(x) = {xの木の根}
        if (par[x] == x) return x;
        return par[x] = root(par[x]);
    }
 
    void unite(int x, int y) { // xとyの木を併合
        int rx = root(x); //xの根をrx
        int ry = root(y); //yの根をry
        if (rx == ry) return; //xとyの根が同じ(=同じ木にある)時はそのまま
        par[rx] = ry; //xとyの根が同じでない(=同じ木にない)時：xの根rxをyの根ryにつける
    }
 
    bool same(int x, int y) { // 2つのデータx, yが属する木が同じならtrueを返す
        int rx = root(x);
        int ry = root(y);
        return rx == ry;
    }
};
 
 
 
int main () {
 
    int n,m,k;cin>>n>>m>>k;
 
    UnionFind tree(n);
 
    priority_queue<Pp> pq;
 
    rep(i,m) {
        int a,b,c;cin>>a>>b>>c;
        a--;b--;
        pq.push(Pp((-1)*c,P(a,b)));
    }
 
    int cnt = 0;
 
    ll ans=0;
 
    while (cnt<n-k) {
        Pp d = pq.top();pq.pop();
        int a = d.second.first;
        int b = d.second.second;
        int c = d.first*(-1);
        if (tree.root(a)!=tree.root(b)) {
            cnt++;
            tree.unite(a,b);
            ans+=c;
        }
    }
 
    cout<<ans<<endl;
 
 
 
}