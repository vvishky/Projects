#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <mpi.h>

using namespace std;

// Function to print a vector
void printVector(const vector<double> &v, const string &name)
{
    cout << name << ": [";
    for (size_t i = 0; i < v.size(); ++i)
    {
        cout << fixed << setprecision(6) << v[i];
        if (i < v.size() - 1)
            cout << ", ";
    }
    cout << "]" << endl;
}

// Function to print a matrix
void printMatrix(const vector<vector<double>> &mat, const string &name)
{
    cout << name << ":\n";
    for (const auto &row : mat)
    {
        for (double val : row)
        {
            cout << fixed << setprecision(6) << val << " ";
        }
        cout << endl;
    }
}

int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    int n, rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get current process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get number of processes

    double start_time, end_time, time_taken;

    vector<vector<double>> L = {
        {0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 1, 1},
        {1, 0, 0, 0, 1, 0},
        {1, 1, 0, 0, 0, 1},
        {0, 0, 1, 0, 1, 0}};

    // Get the number of rows
    n = L.size(); // Number of webpages

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<double> P_flat(n * n); // 1D flattened matrix P
    vector<double> q(n);          // Resultant vector q

    // Compute number of outgoing links (n_j)
    vector<double> n_outgoing(n, 0.0);
    for (int j = 0; j < n; ++j)
    {
        for (int i = 0; i < n; ++i)
        {
            n_outgoing[j] += L[i][j];
        }
    }

    // Create Q matrix
    vector<vector<double>> Q(n, vector<double>(n, 0.0));
    for (int j = 0; j < n; ++j)
    {
        if (n_outgoing[j] > 0)
        {
            for (int i = 0; i < n; ++i)
            {
                Q[i][j] = L[i][j] / n_outgoing[j];
            }
        }
    }

    // Check for dead-end pages (columns all zero)
    vector<double> d(n, 0.0);
    for (int j = 0; j < n; ++j)
    {
        bool allZero = true;
        for (int i = 0; i < n; ++i)
        {
            if (Q[i][j] != 0)
            {
                allZero = false;
                break;
            }
        }
        d[j] = allZero ? 1.0 : 0.0;
    }

    // Create e (vector of ones)
    vector<double> e(n, 1.0);

    // Compute correction for dead-end pages (e*d^t)/n
    vector<vector<double>> deadEnd(n, vector<double>(n, 0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            deadEnd[i][j] = (d[j] * e[i]) / n;
        }
    }

    // Create P (stochastic matrix) with correction for dead-end pages
    vector<vector<double>> P(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            P[i][j] = Q[i][j] + deadEnd[i][j];
        }
    }

    // Flatten the 2D matrix P into a 1D vector for MPI
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            P_flat[i * n + j] = P[i][j];
        }
    }

    // Define row-wise scatter/gather parameters
    vector<int> sendcounts(size), displs(size);
    vector<int> recvcounts(size), recvdispls(size);

    // Row distribution for scatter/gather
    int rows = n, base = rows / size, rem = rows % size;
    int offset = 0;
    for (int i = 0; i < size; ++i)
    {
        int rows_i = base + (i < rem ? 1 : 0);
        sendcounts[i] = rows_i * n;
        displs[i] = offset;
        recvcounts[i] = rows_i;
        recvdispls[i] = (i == 0) ? 0 : recvdispls[i - 1] + recvcounts[i - 1];
        offset += sendcounts[i];
    }

    int local_rows = recvcounts[rank];
    vector<double> local_P(local_rows * n); // Local part of flattened P
    vector<double> local_q(local_rows);     // Local result of matrix-vector multiplication

    // Scatter matrix P among processes
    MPI_Scatterv(P_flat.data(), sendcounts.data(), displs.data(), MPI_DOUBLE,
                 local_P.data(), sendcounts[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    // Power Method with Rayleigh Quotient
    vector<double> r(n, 1.0 / n); // Initial guess (uniform)

    int max_iter = 30;
    double final_rayleigh = 0.0;

    // Broadcast initial r to all processes
    MPI_Bcast(r.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    start_time = MPI_Wtime();
    // Power iteration loop Parallelization with MPI
    for (int iter = 0; iter < max_iter; ++iter)
    {

        // Matrix-vector multiplication (local part)
        for (int i = 0; i < local_rows; ++i)
        {
            local_q[i] = 0;
            for (int j = 0; j < n; ++j)
            {
                local_q[i] += local_P[i * n + j] * r[j];
            }
        }


        // Compute L1 norm of q for normalization
        double local_norm = 0.0;
        for (int i = 0; i < local_rows; ++i)
        {
            local_norm += fabs(local_q[i]);
        }

        // Reduce to get global norm
        double global_norm = 0.0;
        MPI_Allreduce(&local_norm, &global_norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        // Normalize local q to get local r_new
        vector<double> local_r_new(local_rows);
        for (int i = 0; i < local_rows; ++i)
            local_r_new[i] = local_q[i] / global_norm;

        // Gather r_new to rank 0
        vector<double> r_new;
        if (rank == 0)
        {
            r_new.resize(n);
        }
        
        MPI_Gatherv(local_r_new.data(), local_rows, MPI_DOUBLE,
                    r_new.data(), recvcounts.data(), recvdispls.data(), MPI_DOUBLE,
                    0, MPI_COMM_WORLD);

        // Compute local Rayleigh components
        double local_dot_product = 0.0;
        double local_l2norm_sq = 0.0;
        for (int i = 0; i < local_rows; ++i)
        {
            local_dot_product += local_r_new[i] * local_q[i];
            local_l2norm_sq += local_r_new[i] * local_r_new[i];
        }

        // Reduce to get global Rayleigh quotient
        double global_dot_product = 0.0;
        double global_l2norm_sq = 0.0;
        MPI_Reduce(&local_dot_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&local_l2norm_sq, &global_l2norm_sq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // Only rank 0 updates r and Rayleigh
        if (rank == 0)
        {

            final_rayleigh = global_dot_product / global_l2norm_sq;
            r = r_new;
        }

        // Broadcast updated r to all processes for next iteration
        MPI_Bcast(r.data(), n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // Final results
    if (rank == 0)
    {

        end_time = MPI_Wtime();
        time_taken = end_time - start_time;

        printVector(r, "Final PageRank vector (r)");

        cout << "Final Rayleigh Quotient (approx largest eigenvalue): " << fixed << setprecision(6) << final_rayleigh << endl;
        //cout << "Total computation time with " << size << " processors: " << time_taken << " seconds." << endl;

        printMatrix(P,"Stochastic matrix P");  

        // Verify second condition for Page rank
        cout << "\nVerifying Summation condition:\n";
        double sum_r = 0.0;
        for (double val : r)
            sum_r += val;
        cout << "Sum of Final PageRank vector (r): " << fixed << setprecision(3) << sum_r << endl;
    }

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
