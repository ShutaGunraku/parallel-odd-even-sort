#include <iostream>
#include <algorithm>
#include <chrono>

using namespace std;

// function prototype
int* list1_sort(int A[], int n);
int* list2_sort(int A[], int n);
int* list3_sort(int A[], int n);
int* list4_sort(int A[], int n);
int* list4_sort_parallel(int A[], int n);

int main() {
    int A[] = {5, 2, 4, 6, 1, 7, 3, 5, 2, 4, 6, 1, 7, 9};
    int n = sizeof(A) / sizeof(A[0]);

    // Listing 1
    // int A1[n];
    // std::copy(A, A + n, A1);
    // auto start1 = std::chrono::high_resolution_clock::now();
    // int* A1_sorted = list1_sort(A1, n);
    // auto end1 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> duration1 = end1 - start1;

    // Print the sorted arrays
    // cout << "Sorted arrays:" << endl;
    // // listing 1
    // cout << "A1_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A1_sorted[i] << " ";
    // // }
    // cout << endl;
    // cout << "Duration: " << duration1.count() << " seconds" << endl;


    // Listing 2
    // int A2[n];
    // std::copy(A, A + n, A2);
    // auto start2 = std::chrono::high_resolution_clock::now();
    // int* A2_sorted = list2_sort(A2, n);
    // auto end2 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> duration2 = end2 - start2;

    // // listing 2
    // cout << endl;
    // cout << "A2_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A2_sorted[i] << " ";
    // }
    // cout << endl;
    // cout << "Duration: " << duration2.count() << " seconds" << endl;


    // Listing 3
    // int A3[n];
    // std::copy(A, A + n, A3);
    // auto start3 = std::chrono::high_resolution_clock::now();
    // int* A3_sorted = list3_sort(A3, n);
    // auto end3 = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::milli> duration3 = end3 - start3;

    // // listing 3
    // cout << endl;
    // cout << "A3_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A3_sorted[i] << " ";
    // }
    // cout << endl;
    // cout << "Duration: " << duration3.count() << " seconds" << endl;


    // Listing 4
    int A4[n];
    std::copy(A, A + n, A4);
    auto start4 = std::chrono::high_resolution_clock::now();
    int* A4_sorted = list4_sort(A4, n);
    auto end4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration4 = end4 - start4;

    // listing 4
    cout << endl;
    cout << "A4_sorted: ";
    for (int i = 0; i < n; i++) {
        cout << A4_sorted[i] << " ";
    }
    cout << endl;
    cout << "Duration: " << duration4.count() << " seconds" << endl;

    // Listing 4 parallel
    int A4_parallel[n];
    std::copy(A, A + n, A4_parallel);
    auto start4_parallel = std::chrono::high_resolution_clock::now();
    int* A4_sorted_parallel = list4_sort_parallel(A4_parallel, n);
    auto end4_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration4_parallel = end4_parallel - start4_parallel;

    cout << endl;
    cout << "A4_sorted_parallel: ";
    for (int i = 0; i < n; i++) {
        cout << A4_sorted_parallel[i] << " ";
    }
    cout << endl;
    cout << "Duration: " << duration4_parallel.count() << " seconds" << endl;

    return 0;
}

int* list1_sort(int A[], int n)
{
    // Parallel Odd-Even Sort
    for (int p = 1; p < n; p += p) 
        for (int k = p; k > 0; k /= 2) 
            for (int j = k % p; j + k < n; j += 2 * k) 
                for (int i = 0; i < n-j-k; i++)
                    if ((j+i) / (p+p) == (j+i+k)/(p+p))
                        if (A[j+i] > A[j+i+k])
                            swap(A[j+i], A[j+i+k]);

    return A;
}

int* list2_sort(int A[], int n)
{
    // Parallel Odd-Even Sort
    for(int p = 1; p < n; p *= 2) 
        for(int k = p; k > 0; k /= 2) 
            for(int j = k % p; j + k < 2*p; j += 2*k) 
                for(int i = 0; i < k; i++) 
                    for(int m = i + j; m < n - k; m += 2*p) 
                        if(A[m] > A[m+k]) 
                            swap(A[m], A[m+k]);

    return A;
}

int* list3_sort(int A[], int n)
{
    for (int p = 1; p < n; p *= 2) 
        for (int k = p; k > 0; k /= 2) 
            for (int j = k % p; j + k < n; j += 2*k) 
                for (int i = std::min(k, n-j-k); i--;) 
                    if ((j+i)/(2*p) == (j+i+k)/(2*p)) 
                        if (A[j+i] > A[j+i+k]) 
                            std::swap(A[j+i], A[j+i+k]);

    return A;
}

int* list4_sort(int A[], int n)
{
    for(int p = 1; p < n; p *= 2)
        for(int k = p; k > 0; k /= 2)
            for(int j = k & (p - 1); j + k < n; j += 2*k)
                if((j | (2*p - 1)) == ((j+k) | (2*p - 1)))
                    for(int i = std::min(k, n-j-k); i--;)
                        if(A[j+i] > A[j+i+k])
                            std::swap(A[j+i], A[j+i+k]);

    return A;
}

int* list4_sort_parallel(int A[], int n)
{
    /*
    "pragma omp parallel for" arallelise the outermost loop across multiple threads.
    "default(none)" and "shared(A,n)" specifies that the variable
    A and n are shared between the threads.
    */
    #pragma omp parallel for default(none) shared(A,n)
    for(int p = 1; p < n; p *= 2)
    {
        /*
        "pragma omp for" statements parallelise the loops over k, j and i.
        */
        #pragma omp for
        for(int k = p; k > 0; k /= 2)
        {
            #pragma omp for
            for(int j = k & (p - 1); j + k < n; j += 2*k)
            {
                if((j | (2*p - 1)) == ((j+k) | (2*p - 1)))
                {
                    #pragma omp for
                    for(int i = std::min(k, n-j-k); i--;)
                    {
                        if(A[j+i] > A[j+i+k])
                        {
                            std::swap(A[j+i], A[j+i+k]);
                        }
                    }
                }
            }
        }
    }

    return A;
}