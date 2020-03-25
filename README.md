# competitive-programing

#### Union Find

最初n個の無関係な点を関係性ごとにつないでいける

友達の友達は友達

[https://atcoder.jp/contests/arc065/tasks/arc065_b](https://atcoder.jp/contests/arc065/tasks/arc065_b)

UnionFind 2個用意して、道路と鉄道の情報を管理して、

n個の点それぞれの、道路のrootと鉄道のrootを出して、mapでその組み合わせの

出現回数を数える





#### BFS

幅優先探索

全ての点から点が等距離の場合は、最短距離を求めることができる

https://atcoder.jp/contests/arc005/tasks/arc005_3

'.'は通れて、'#'は通れないが、二回までなら'#'を通っても良い

'.'の時、cost=0,'#'の時、cost=1として、それぞれのますに行くためのcostの最小を更新していく。

cost(gx,gy)<=2ならok

<details><summary>コード</summary><div>
    
```
#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
const long long INF = 1LL<<60;
#define rev(s) (string((s).rbegin(), (s).rend()))

int h,w,sy,sx,gy,gx;;
const int inf = 300000;
vector<string> field(510);
int dx[] = {1,0,-1,0};
int dy[] = {0,1,0,-1};
vector<vector<int>> cost(510,vector<int>(510,inf));

void bfs() {
    int cs;
    queue<P> q;
    q.push(P(sy,sx));
    cost[sy][sx]=0;
    while (q.size()) {
        P pa = q.front();q.pop();
        rep(i,4) {
            int ny = pa.first+dy[i];
            int nx = pa.second+dx[i];
            if (ny>=0&&ny<h&&nx>=0&&nx<w) {
                if (field[ny][nx]=='#') {
                    cs=1;
                }
                else {
                    cs=0;
                }
                if (cost[ny][nx]>cost[pa.first][pa.second]+cs) {
                    cost[ny][nx]=cost[pa.first][pa.second]+cs;
                    q.push(P(ny,nx));
                }
            }
        }
    }
}


int main() {
    cin>>h>>w;
    rep(i,h) {
        cin>>field[i];
    }
    rep(i,h) {
        rep(j,w) {
            if (field[i][j]=='s') {
                sy=i;sx=j;
            }if (field[i][j]=='g') {
                gy=i;gx=j;
            }
        }
    }
    bfs();
    if (cost[gy][gx]<=2) {
        cout<<"YES"<<endl;
    }
    else {
        cout<<"NO"<<endl;
    }

}

```

</div></details>



https://atcoder.jp/contests/agc043/tasks/agc043_a

上の問題の応用

