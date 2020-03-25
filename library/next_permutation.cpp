#include <iostream>
#include <algorithm>
using namespace std;

int N, A[12][12], B[12], perm[12], ans = 2000000000;

int main(){
    int array[]={1,2,3};
    do{
        for(int i=0; i<3; i++){
            cout<<array[i];
            if(i!=2)cout<<" ";
        }
        cout<<endl;
    }while(next_permutation(array,array+3));
    return 0;
}