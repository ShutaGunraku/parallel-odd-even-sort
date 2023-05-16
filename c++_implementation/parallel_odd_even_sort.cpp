/**
 * @file parallel_odd_even_sort.cpp
 * @brief This program generates arrays with random integers, sorts them using different sorting functions, and exports the execution times to a CSV file.
 * @version 1.0
 * @date 14th May 2023
 * @author Shuta Gunraku
 */

#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

// Function prototype
void exportToCSV(const std::string& filename, const std::vector<long double>& data, const std::vector<long>& sizes);

long double executeListing(int A[], int n, int* (*sortFunc)(int[], int), const std::string& sortFuncName);

int* sortListing1(int A[], int n);
int* sortListing1Parallel(int A[], int n);

int* sortListing2(int A[], int n);
int* sortListing2Parallel(int A[], int n);

int* sortListing3(int A[], int n);
int* sortListing3Parallel(int A[], int n);

int* sortListing4(int A[], int n);
int* sortListing4Parallel(int A[], int n);

/**
 * Main function of the program
 * @return 0 if the program is successful
*/
int main() 
{
    // Generate random numbers
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(1, 100); // Configure the range of the random numbers

    // Configure the size of arrays to be sorted
    std::vector<long> sizes;
    sizes.push_back(8);
    sizes.push_back(9);


    std::vector<long double> executionTimes;

    // Generate arrays of random numbers and sort them
    for (int zeroCount: sizes) {
        int n = pow(10, zeroCount);
        int* A = new int[n];
        for (int i = 0; i < n; i++) {
            A[i] = distribution(generator);
        }

        cout << "Sorting an array of length n = pow(10," << zeroCount << ")" << endl;

        // Function pointers
        auto sortFunc1 = sortListing1;
        auto sortFunc1Parallel = sortListing1Parallel;
        auto sortFunc2 = sortListing2;
        auto sortFunc2Parallel = sortListing2Parallel;
        auto sortFunc3 = sortListing3; 
        auto sortFunc3Parallel = sortListing3Parallel;
        auto sortFunc4 = sortListing4;
        auto sortFunc4Parallel = sortListing4Parallel;
        
        // Execute the sorting functions and measure execution time
        // executionTimes.push_back(executeListing(A, n, sortFunc1, "Listing1"));
        // executionTimes.push_back(executeListing(A, n, sortFunc1Parallel, "Listing1Parallel"));
        executionTimes.push_back(executeListing(A, n, sortFunc2, "Listing2"));
        executionTimes.push_back(executeListing(A, n, sortFunc2Parallel, "Listing2Parallel"));
        // executionTimes.push_back(executeListing(A, n, sortFunc3, "Listing3"));
        // executionTimes.push_back(executeListing(A, n, sortFunc3Parallel, "Listing3Parallel"));
        // executionTimes.push_back(executeListing(A, n, sortFunc4, "Listing4"));
        // executionTimes.push_back(executeListing(A, n, sortFunc4Parallel, "Listing4Parallel"));

        // Deallocate the dynamic array
        delete[] A;
    }

    // Export the vector to a CSV file
    exportToCSV("output.csv", executionTimes, sizes);

    return 0;
}


/**
 * Exports the execution times to a CSV file.
 * 
 * @param filename the name of the CSV file
 * @param data the execution times
 * @param sizes the sizes of arrays to be sorted
 * @return void
*/
void exportToCSV(const std::string& filename, const std::vector<long double>& data, const std::vector<long>& sizes)
{
    std::ofstream outputFile(filename);
    if (!outputFile.is_open())
    {
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write column headers
    outputFile << "n,Listing1,Listing1Parallel,Listing2,Listing2Parallel,Listing3,Listing3Parallel,Listing4,Listing4Parallel" << std::endl;

    int index = 0;
    for (int n : sizes)
    {
        outputFile << n << ",";
        for (int i = 0; i < 8; i++)
        {
            outputFile << data[index] << ",";
            index++;
        }
        outputFile << std::endl;
    }

    outputFile.close();
}


/**
 * Executes a sorting function and measures the execution time.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @param sortFunc the sorting function
 * @param sortFuncName the name of the sorting function
 * @return the execution time
 * @see https://stackoverflow.com/questions/22387586/measuring-execution-time-of-a-function-in-c
*/
long double executeListing(int A[], int n, int* (*sortFunc)(int[], int), const std::string& sortFuncName)
{
    // Copy the array
    int* ACopied = new int[n];
    std::copy(A, A + n, ACopied);

    // Execute the sorting function and measure execution time
    auto start = std::chrono::high_resolution_clock::now();
    int* ASorted = sortFunc(ACopied, n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    long double executionTime = duration.count() / 1000.0;

    // Print the sorted array and execution time
    cout << sortFuncName << ": ";
    // for (int i = 0; i < n; i++) {
    //     cout << ASorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << executionTime << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] ACopied;

    return executionTime;
}


/**
 * Sorts an array of integers using Listing 1.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
 * @see R. Sedgewick, "Algorithms in C++, 1992," ed: Addison-Wesley.
*/
int* sortListing1(int A[], int n)
{
    for (int p = 1; p < n; p += p) 
        for (int k = p; k > 0; k /= 2) 
            for (int j = k % p; j + k < n; j += 2 * k) 
                for (int i = 0; i < n-j-k; i++)
                    if ((j+i) / (p+p) == (j+i+k)/(p+p))
                        if (A[j+i] > A[j+i+k])
                            swap(A[j+i], A[j+i+k]);

    return A;
}


/**
 * @brief Sorts an array of integers using Listing 1 in parallel.
 * "#pragma omp parallel for" parallelises for loops.
 * "shared()" means that the variables inside it are shared between the threads.
 * "default(none)" means that all variables must be explicitly declared as shared or private.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing1Parallel(int A[], int n)
{

    for (int p = 1; p < n; p *= 2) 
    {
        #pragma omp parallel for shared(A, n, p) default(none)
        for (int k = p; k > 0; k /= 2) 
        {
            #pragma omp parallel for shared(A, n, p, k) default(none)
            for (int j = k % p; j + k < n; j += 2 * k) 
            {
                #pragma omp parallel for shared(A, n, p, k, j) default(none)
                for (int i = 0; i < n-j-k; i++)
                {
                    if ((j+i) / (p+p) == (j+i+k)/(p+p))
                    {
                        if (A[j+i] > A[j+i+k])
                        {
                            // Ensures atomic access to the shared variable
                            #pragma omp critical
                            swap(A[j+i], A[j+i+k]);
                        }
                    }
                }
            }
        }
    }

    return A;
}


/**
 * Sorts an array of integers using Listing 2.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing2(int A[], int n)
{
    for(int p = 1; p < n; p *= 2) 
        for(int k = p; k > 0; k /= 2) 
            for(int j = k % p; j + k < 2*p; j += 2*k) 
                for(int i = 0; i < k; i++) 
                    for(int m = i + j; m < n - k; m += 2*p) 
                        if(A[m] > A[m+k]) 
                            swap(A[m], A[m+k]);

    return A;
}


/**
 * @brief Sorts an array of integers using Listing 2 in parallel.
 * "#pragma omp parallel for" parallelises for loops.
 * "shared()" means that the variables inside it are shared between the threads.
 * "default(none)" means that all variables must be explicitly declared as shared or private.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
// int* sortListing2Parallel(int A[], int n)
// {
//     // no parallelisation part
//     for(int p = 1; p < n; p *= 2) 
//     {
//         #pragma omp parallel for shared(A, n, p) default(none)
//         for(int k = p; k > 0; k /= 2) 
//         {
//             #pragma omp parallel for shared(A, n, p, k) default(none)
//             for(int j = k % p; j + k < 2*p; j += 2*k) 
//             {
//                 #pragma omp parallel for shared(A, n, p, k, j) default(none)
//                 for(int i = 0; i < k; i++) 
//                 {
//                     #pragma omp parallel for shared(A, n, p, k, j, i) default(none)
//                     for(int m = i + j; m < n - k; m += 2*p) 
//                     {
//                         if(A[m] > A[m+k]) 
//                         {
//                             // "#pragma omp critical" ensures that only one thread can access the shared variable at a time.
//                             #pragma omp critical
//                             swap(A[m], A[m+k]);
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     return A;
// }
int* sortListing2Parallel(int A[], int n)
{
    for(int p = 1; p < n; p *= 2) 
    {
        #pragma omp parallel for shared(A, n, p) default(none)
        for(int k = p; k > 0; k /= 2) 
        {
            #pragma omp parallel for shared(A, n, p, k) default(none)
            for(int j = k % p; j + k < 2*p; j += 2*k) 
            {
                #pragma omp parallel for shared(A, n, p, k, j) default(none)
                for(int i = 0; i < k; i++) 
                {
                    #pragma omp parallel for shared(A, n, p, k, j, i) reduction(min : m) default(none)
                    for(int m = i + j; m < n - k; m += 2*p) 
                    {
                        if(A[m] > A[m+k]) 
                        {
                            swap(A[m], A[m+k]);
                        }
                    }
                }
            }
        }
    }
    return A;
}


/**
 * Sorts an array of integers using Listing 3.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing3(int A[], int n)
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


/**
 * @brief Sorts an array of integers using Listing 3 in parallel.
 * "#pragma omp parallel for" parallelises for loops.
 * "shared()" means that the variables inside it are shared between the threads.
 * "default(none)" means that all variables must be explicitly declared as shared or private.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing3Parallel(int A[], int n)
{
    for (int p = 1; p < n; p *= 2) 
    {
        #pragma omp parallel for shared(A, n, p) default(none)
        for (int k = p; k > 0; k /= 2) 
        {
            #pragma omp parallel for shared(A, n, p, k) default(none)
            for (int j = k % p; j + k < n; j += 2*k) 
            {
                #pragma omp parallel for shared(A, n, p, k, j) default(none)
                for (int i = std::min(k, n-j-k); i--;) 
                {
                    if ((j+i)/(2*p) == (j+i+k)/(2*p)) 
                    {
                        if (A[j+i] > A[j+i+k]) 
                        {
                            // Ensures atomic access to the shared variable
                            #pragma omp critical
                            std::swap(A[j+i], A[j+i+k]);
                        }
                    }
                }
            }
        }
    }

    return A;
}


/**
 * Sorts an array of integers using Listing 4.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing4(int A[], int n)
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


/**
 * @brief Sorts an array of integers using Listing 4 in parallel.
 * "#pragma omp parallel for" parallelises for loops.
 * "shared()" means that the variables inside it are shared between the threads.
 * "default(none)" means that all variables must be explicitly declared as shared or private.
 * 
 * @param A the array to be sorted
 * @param n the size of the array
 * @return the sorted array
*/
int* sortListing4Parallel(int A[], int n)
{
    #pragma omp parallel for default(none) shared(A,n)
    for(int p = 1; p < n; p *= 2)
    {
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