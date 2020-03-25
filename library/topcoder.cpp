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

class MultiplicationTable2{
    public:
    int minimalGoodSet(vector <int> table){
        int n = sqrt(table.size());
        int vec[n][n];
        rep(i,n) {
            rep(j,n) {
                vec[i][j]=table[n*i+j];
            }
        }
        int ans=2500;
        rep(i,n) {
            for (int j=0;j<=n-1-i;j++) {
                bool ok=true;
                rep(k,j+1) {
                    rep(l,j+1) {
                        if (i>vec[i+k][i+l] || i+j<vec[i+k][i+l]) {
                            ok=false;
                            break;
                        }
                    }
                    if (!ok) break;
                }
                if (ok) {
                    ans=min(ans,i+j+1);
                }
            }
        }
        return ans;
    }
};

int main() {

    MultiplicationTable2 mt;

    int N;cin>>N;

    vector<int> ary;

    rep(i,N) {
        int p;cin>>p;
        ary.push_back(p);
    }

    int ans = mt.minimalGoodSet(ary);

    cout<<ans<<endl;

}