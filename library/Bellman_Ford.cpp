#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
typedef pair<int,int> P;
const long long INF = 1LL<<60;
#define rev(s) (string((s).rbegin(), (s).rend()))
// cout << fixed << setprecision(10) << ans << endl; 有効桁数指定

int main () {
    int n,m;cin>>n>>m;
    ll p;cin>>p;
    
    int a[5010],b[5010];
    ll c[5010];
    rep(i,m) {
        cin>>a[i]>>b[i]>>c[i];
        a[i]--;
        b[i]--;
        c[i]=-(c[i]-p);
    }

    ll dist[2510];
    rep(i,n) {
        dist[i]=INF;
    }
    dist[0]=0;

    rep(i,n-1) {
        rep(j,m) {
            if (dist[a[j]]==INF) continue;

            if (dist[b[j]]>dist[a[j]]+c[j]) {
                dist[b[j]]=dist[a[j]]+c[j];
            }
        }
    }

    ll ans = dist[n-1];

    bool negative[2510];
    rep(i,n) negative[i]=false;

    rep(i,n) {
        rep(j,m) {
            if (dist[a[j]]==INF) continue;

            if (dist[b[j]]>dist[a[j]]+c[j]) {
                dist[b[j]]=dist[a[j]]+c[j];
                negative[j]=true;
            }
            if (negative[a[j]]) negative[b[j]]==true;
        }
    }

    if (negative[n-1]) cout<<-1<<endl;
    else {
        if (ans>0) {
            cout<<0<<endl;
        }
        else {
            cout<<-ans<<endl;
        }
    }

}