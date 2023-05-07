#include <iostream>
#include <algorithm>
#include <chrono>

using namespace std;

// function prototype
int* list1_sort(int A[], int n);
int* list2_sort(int A[], int n);
int* list3_sort(int A[], int n);

int main() {
    int A[] = {5, 2, 4, 6, 1, 7, 3, 5, 2, 4, 6, 1, 7, 9};
    int n = sizeof(A) / sizeof(A[0]);

    // Listing 1
    int A1[n];
    std::copy(A, A + n, A1);
    auto start1 = std::chrono::high_resolution_clock::now();
    int* A1_sorted = list1_sort(A1, n);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration1 = end1 - start1;

    // Listing 2
    int A2[n];
    std::copy(A, A + n, A2);
    auto start2 = std::chrono::high_resolution_clock::now();
    int* A2_sorted = list2_sort(A2, n);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration2 = end2 - start2;

    // Listing 3
    int A3[n];
    std::copy(A, A + n, A3);
    auto start3 = std::chrono::high_resolution_clock::now();
    int* A3_sorted = list3_sort(A3, n);
    auto end3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration3 = end3 - start3;

    // Print the sorted arrays
    cout << "Sorted arrays:" << endl;
    // listing 1
    cout << "A1_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A1_sorted[i] << " ";
    }
    cout << endl;
    cout << "Duration: " << duration1.count() << " seconds" << endl;

    // listing 2
    cout << endl;
    cout << "A2_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A2_sorted[i] << " ";
    }
    cout << endl;
    cout << "Duration: " << duration2.count() << " seconds" << endl;

    // listing 3
    cout << endl;
    cout << "A3_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A3_sorted[i] << " ";
    }
    cout << endl;
    cout << "Duration: " << duration3.count() << " seconds" << endl;

    return 0;
}

int* list1_sort(int A[], int n)
{
    // Parallel Odd-Even Sort
    for (int p = 1; p < n; p += p) {
        for (int k = p; k > 0; k /= 2) {
            for (int j = k % p; j + k < n; j += 2 * k) {
                for (int i = 0; i < n-j-k; i++){
                    if ((j+i) / (p+p) == (j+i+k)/(p+p)){
                        if (A[j+i] > A[j+i+k]){
                            swap(A[j+i], A[j+i+k]);
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

int* list3_sort(int A[], int n)
{
    for (int p = 1; p < n; p *= 2) {
        for (int k = p; k > 0; k /= 2) {
            for (int j = k % p; j + k < n; j += 2*k) {
                for (int i = std::min(k, n-j-k); i--;) {
                    if ((j+i)/(2*p) == (j+i+k)/(2*p)) {
                        if (A[j+i] > A[j+i+k]) {
                            std::swap(A[j+i], A[j+i+k]);
                        }
                    }   
                }
            }
        }
    }

    return A;
}