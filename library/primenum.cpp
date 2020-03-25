#include <iostream>
#include <vector>
using namespace std;

int main() {
    // エラトステネスのふるい
    int MAX = 101010;
    vector<int> is_prime(MAX, 1);
    is_prime[0] = 0, is_prime[1] = 0;
    for (int i = 2; i < MAX; ++i) {
        if (!is_prime[i]) continue;
        for (int j = i*2; j < MAX; j += i) is_prime[j] = 0;
    }
}