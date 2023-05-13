#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>


using namespace std;


// function prototype
void exportToCSV(const std::string& filename, const std::vector<long double>& data);

long double executeListing1(int A[], int n);
long double executeListing1Parallel(int A[], int n);
int* sortListing1(int A[], int n);
int* sortListing1Parallel(int A[], int n);

long double executeListing2(int A[], int n);
long double executeListing2Parallel(int A[], int n);
int* sortListing2(int A[], int n);
int* sortListing2Parallel(int A[], int n);

long double executeListing3(int A[], int n);
long double executeListing3Parallel(int A[], int n);
int* sortListing3(int A[], int n);
int* sortListing3Parallel(int A[], int n);

long double executeListing4(int A[], int n);
long double executeListing4Parallel(int A[], int n);
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

    long double A1_time = executeListing1(A, n);
    long double A1_parallel_time = executeListing1Parallel(A, n);

    long double A2_time = executeListing2(A, n);
    long double A2_parallel_time = executeListing2Parallel(A, n);

    long double A3_time = executeListing3(A, n);
    long double A3_sorted_time = executeListing3Parallel(A, n);

    long double A4_sorted = executeListing4(A, n);
    long double A4_sorted_time = executeListing4Parallel(A, n);

    // Deallocate dynamic arrays
    delete[] A;

    // Store the long double values in a vector
    std::vector<long double> values;
    values.push_back(A1_time);
    values.push_back(A1_parallel_time);
    values.push_back(A2_time);
    values.push_back(A2_parallel_time);
    values.push_back(A3_time);
    values.push_back(A3_sorted_time);
    values.push_back(A4_sorted);
    values.push_back(A4_sorted_time);

    // Export the vector to a CSV file
    exportToCSV("output.csv", values);

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

// Listing 1
long double executeListing1(int A[], int n)
{
    int* A1 = new int[n];
    std::copy(A, A + n, A1);
    auto start1 = std::chrono::high_resolution_clock::now();
    int* A1_sorted = sortListing1(A1, n);
    auto end1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration1 = end1 - start1;

    long double execution_time = duration1.count() / 1000.0;
    // Print the sorted arrays
    cout << "Sorted arrays:" << endl;
    // listing 1
    cout << "A1_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A1_sorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A1;

    return execution_time;
}

// Listing 1 parallel
long double executeListing1Parallel(int A[], int n)
{
    int* A1_parallel = new int[n];
    std::copy(A, A + n, A1_parallel);
    auto start1_parallel = std::chrono::high_resolution_clock::now();
    int* A1_sorted_parallel = sortListing1Parallel(A1_parallel, n);
    auto end1_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration1_parallel = end1_parallel - start1_parallel;

    long double execution_time = duration1_parallel.count() / 1000.0;
    cout << endl;
    cout << "A1_sorted_parallel: ";
    // for (int i=0; i<n; i++) {
    //     cout << A1_sorted_parallel[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A1_parallel;

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


// Listing 2
long double executeListing2(int A[], int n)
{
    int* A2 = new int[n];
    std::copy(A, A + n, A2);
    auto start2 = std::chrono::high_resolution_clock::now();
    int* A2_sorted = sortListing2(A2, n);
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration2 = end2 - start2;

    long double execution_time = duration2.count() / 1000.0;
    cout << endl;
    cout << "A2_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A2_sorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A2;

    return execution_time;
}


// listing 2 parallel
long double executeListing2Parallel(int A[], int n)
{
    int* A2_parallel = new int[n];
    std::copy(A, A + n, A2_parallel);
    auto start2_parallel = std::chrono::high_resolution_clock::now();
    int* A2_sorted_parallel = sortListing2Parallel(A2_parallel, n);
    auto end2_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration2_parallel = end2_parallel - start2_parallel;

    long double execution_time = duration2_parallel.count() / 1000.0;
    cout << endl;
    cout << "A2_sorted_parallel: ";
    // for (int i=0; i<n; i++) {
    //     cout << A2_sorted_parallel[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A2_parallel;

    return execution_time;
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

// Listing 3
long double executeListing3(int A[], int n)
{
    int* A3 = new int[n];
    std::copy(A, A + n, A3);
    auto start3 = std::chrono::high_resolution_clock::now();
    int* A3_sorted = sortListing3(A3, n);
    auto end3 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration3 = end3 - start3;

    long double execution_time = duration3.count() / 1000.0;
    cout << endl;
    cout << "A3_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A3_sorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A3;

    return execution_time;
}

long double executeListing3Parallel(int A[], int n)
{
    // Listing 3 parallel
    int* A3_parallel = new int[n];
    std::copy(A, A + n, A3_parallel);
    auto start3_parallel = std::chrono::high_resolution_clock::now();
    int* A3_sorted_parallel = sortListing3Parallel(A3_parallel, n);
    auto end3_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration3_parallel = end3_parallel - start3_parallel;

    long double execution_time = duration3_parallel.count() / 1000.0;
    cout << endl;
    cout << "A3_sorted_parallel: ";
    // for (int i=0; i<n; i++) {
    //     cout << A3_sorted_parallel[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A3_parallel;

    return execution_time;
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


// Listing 4
long double executeListing4(int A[], int n)
{
    int* A4 = new int[n];
    std::copy(A, A + n, A4);
    auto start4 = std::chrono::high_resolution_clock::now();
    int* A4_sorted = sortListing4(A4, n);
    auto end4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration4 = end4 - start4;

    long double execution_time = duration4.count() / 1000.0;
    cout << endl;
    cout << "A4_sorted: ";
    // for (int i = 0; i < n; i++) {
    //     cout << A4_sorted[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A4;

    return execution_time;
}

long double executeListing4Parallel(int A[], int n)
{
    // Listing 4 parallel
    int* A4_parallel = new int[n];
    std::copy(A, A + n, A4_parallel);
    auto start4_parallel = std::chrono::high_resolution_clock::now();
    int* A4_sorted_parallel = sortListing4Parallel(A4_parallel, n);
    auto end4_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration4_parallel = end4_parallel - start4_parallel;

    long double execution_time = duration4_parallel.count() / 1000.0;
    cout << endl;
    cout << "A4_sorted_parallel: ";
    // for (int i=0; i<n; i++) {
    //     cout << A4_sorted_parallel[i] << " ";
    // }
    cout << endl;
    cout << "Duration: " << execution_time << " seconds" << endl;

    // Deallocate dynamic arrays
    delete[] A4_parallel;

    return execution_time;
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
