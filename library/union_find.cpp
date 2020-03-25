#include <bits/stdc++.h>
#define rep(i,n)for(int i=0;i<(n);i++)
using namespace std;
typedef long long ll;
const long long INF = 1LL<<60;
#define rev(s) (string((s).rbegin(), (s).rend()))
typedef pair<ll,int> p;

vector<int> par(1000001);
vector<int> r(1000001);

void init(int n) {
    rep(i,n) {
        par[i]=i;
        r[i]=0;
    }
}

int root(int x) {
    if (par[x]==x) {
        return x;
    }
    else {
        return par[x]=root(par[x]);
    }
}

void unit(int x,int y) {
    x = root(x);
    y = root(y);
    if (x==y) return;

    else {
        if (r[x]<r[y]) {
            par[x]=y;
        }
        else {
            par[y]=x;
            if (root(x)==root(y)) r[x]++;
        }
    }
}

bool same (int x,int y) {
    return root(x)==root(y);
}

int main () {

    int n,q;cin>>n>>q;
    
    init(n);

    rep(i,q) {
        int p,a,b;cin>>p>>a>>b;

        if (p==0) {
            unit(a,b);
        }
        else {
            if (same(a,b)) {
                cout<<"Yes"<<endl;
            }
            else {
                cout<<"No"<<endl;
            }
        }
    }

}

