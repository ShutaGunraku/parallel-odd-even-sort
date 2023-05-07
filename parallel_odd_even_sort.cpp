#include <iostream>
#include <algorithm>

using namespace std;

// function prototype
int list1_sort(int A[], int n);

int main() {
    int A[] = {5, 2, 4, 6, 1, 7, 3, 5, 2, 4, 6, 1, 7, 9};
    int n = sizeof(A) / sizeof(A[0]);

    int A1[n];
    std::copy(A, A + n, A1);

    list1_sort(A1, n);

    return 0;
}

int list1_sort(int A[], int n)
{
        // Parallel Odd-Even Sort
    for (int p = 1; p < n; p += p) {
        for (int k = p; k > 0; k /= 2) {
            for (int j = k & (p - 1); j + k < n; j += 2 * k) {
                if ((j | (2 * p - 1)) == ((j + k) | (2 * p - 1))) {
                    for (int i = min(k, n - j - k) - 1; i >= 0; i--) {
                        if (A[j + i] > A[j + i + k]) {
                            swap(A[j + i], A[j + i + k]);
                        }
                    }
                }
            }
        }
    }

    // Print the sorted array
    for (int i = 0; i < n; i++) {
        cout << A[i] << " ";
    }

    return 0;
}
