#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>


using namespace std;

// function prototype
void exportToCSV(const std::string& filename, const std::vector<long double>& data);

long double executeListing(int A[], int n, int* (*sortFunc)(int[], int), const std::string& sortFuncName);

int* sortListing1(int A[], int n);
int* sortListing1Parallel(int A[], int n);

int* sortListing2(int A[], int n);
int* sortListing2Parallel(int A[], int n);

int* sortListing3(int A[], int n);
int* sortListing3Parallel(int A[], int n);

int* sortListing4(int A[], int n);
int* sortListing4Parallel(int A[], int n);


int main() {
    // Generate random numbers
    std::random_device rd;
    std::mt19937 generator(rd());
    std::uniform_int_distribution<int> distribution(1, 100);  // Adjust the range as needed
    // 200000000
    int n = 1000;  // Adjust the size as needed
    int* A = new int[n];
    for (int i = 0; i < n; i++) {
        A[i] = distribution(generator);
    }

    std::vector<long double> executionTimes;

    // Define the sorting functions
    auto sortFunc1 = sortListing1;
    auto sortFunc1Parallel = sortListing1Parallel;
    auto sortFunc2 = sortListing2;
    auto sortFunc2Parallel = sortListing2Parallel;
    auto sortFunc3 = sortListing3;
    auto sortFunc3Parallel = sortListing3Parallel;
    auto sortFunc4 = sortListing4;
    auto sortFunc4Parallel = sortListing4Parallel;

    cout << "Sorted arrays:" << endl;

    // Execute the sorting functions and measure execution time
    executionTimes.push_back(executeListing(A, n, sortFunc1, "Listing1"));
    executionTimes.push_back(executeListing(A, n, sortFunc1Parallel, "Listing1Parallel"));
    executionTimes.push_back(executeListing(A, n, sortFunc2, "Listing2"));
    executionTimes.push_back(executeListing(A, n, sortFunc2Parallel, "Listing2Parallel"));
    executionTimes.push_back(executeListing(A, n, sortFunc3, "Listing3"));
    executionTimes.push_back(executeListing(A, n, sortFunc3Parallel, "Listing3Parallel"));
    executionTimes.push_back(executeListing(A, n, sortFunc4, "Listing4"));
    executionTimes.push_back(executeListing(A, n, sortFunc4Parallel, "Listing4Parallel"));

    // Deallocate dynamic arrays
    delete[] A;

    // Export the vector to a CSV file
    exportToCSV("output.csv", executionTimes);

    return 0;
}

void exportToCSV(const std::string& filename, const std::vector<long double>& data)
{
    std::ofstream outputFile(filename);
    if (!outputFile.is_open())
    {
        std::cout << "Error opening file: " << filename << std::endl;
        return;
    }

    for (const auto& value : data)
    {
        outputFile << value << ",";
    }

    outputFile.close();
}


long double executeListing(int A[], int n, int* (*sortFunc)(int[], int), const std::string& sortFuncName)
{
    int* A_copy = new int[n];
    std::copy(A, A + n, A_copy);

    auto start = std::chrono::high_resolution_clock::now();
    int* ASorted = sortFunc(A_copy, n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;

    long double execution_time = duration.count() / 1000.0;

    // Print the sorted array and execution time
    cout << sortFuncName << ": ";
    // for (int i = 0; i < n; i++) {
    //     cout << ASorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A_copy;

    return execution_time;
}


int* sortListing1(int A[], int n)
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

int* sortListing1Parallel(int A[], int n)
{

    for (int p = 1; p < n; p *= 2) 
    {
        /* "#pragram omp parallel for" is used for parallelising the outermost loop.
            shared(A, n, p) means that the variables A, n and p are shared between the threads.
            default(none) means that all variables must be explicitly declared as shared or private.
        */
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
                            // "#pragma omp critical" is used to ensure that only one thread can access the shared variable at a time.
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


int* sortListing2(int A[], int n)
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

int* sortListing2Parallel(int A[], int n)
{
    // no parallelisation part
    for(int p = 1; p < n; p *= 2) 
    {
        /*
        "#pragma omp parallel for" is used to parallelise the loop.
        shared(A, n, p) means that the variables A, n and p are shared between the threads.
        default(none) means that all variables must be explicitly specified as shared or private.
        */
        #pragma omp parallel for shared(A, n, p) default(none)
        for(int k = p; k > 0; k /= 2) 
        {
            #pragma omp parallel for shared(A, n, p, k) default(none)
            for(int j = k % p; j + k < 2*p; j += 2*k) 
            {
                #pragma omp parallel for shared(A, n, p, k, j) default(none)
                for(int i = 0; i < k; i++) 
                {
                    #pragma omp parallel for shared(A, n, p, k, j, i) default(none)
                    for(int m = i + j; m < n - k; m += 2*p) 
                    {
                        if(A[m] > A[m+k]) 
                        {
                            /*
                            critical means that the following block of code is executed by only one thread at a time.
                            */
                            #pragma omp critical
                            swap(A[m], A[m+k]);
                        }
                    }
                }
            }
        }
    }
    return A;
}


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


int* sortListing4Parallel(int A[], int n)
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
