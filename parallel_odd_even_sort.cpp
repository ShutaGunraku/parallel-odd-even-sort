#include <iostream>
#include <algorithm>

using namespace std;

// function prototype
int* list1_sort(int A[], int n);
int* list2_sort(int A[], int n);

int main() {
    int A[] = {5, 2, 4, 6, 1, 7, 3, 5, 2, 4, 6, 1, 7, 9};
    int n = sizeof(A) / sizeof(A[0]);

    // Listing 1
    int A1[n];
    std::copy(A, A + n, A1);
    int* A1_sorted = list1_sort(A1, n);

    // Listing 2
    int A2[n];
    std::copy(A, A + n, A2);
    int* A2_sorted = list2_sort(A2, n);

    // Print the sorted arrays
    cout << "Sorted arrays:" << endl;
    cout << "A1_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A1_sorted[i] << " ";
    }
    cout << endl;
    cout << "A2_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A2_sorted[i] << " ";
    }


    return 0;
}

int* list1_sort(int A[], int n)
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

    return A;
}

int* list2_sort(int A[], int n)
{
    // Parallel Odd-Even Sort
    for(int p = 1; p < n; p *= 2) {
        for(int k = p; k > 0; k /= 2) {
            for(int j = k % p; j + k < 2*p; j += 2*k) {
                for(int i = 0; i < k; i++) {
                    for(int m = i + j; m < n - k; m += 2*p) {
                        if(A[m] > A[m+k]) {
                            swap(A[m], A[m+k]);
                        }
                    }
                }
            }  
        }
    }

    return A;
}