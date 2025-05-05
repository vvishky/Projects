# MPI Google PageRank - HPC Programming Project

This repository contains the implementation of a parallelized Google PageRank algorithm using MPI, as part of the *High Performance Computing and Optimization* course at Technische Universität Bergakademie Freiberg.

## 📌 Project Description

The goal of this project is to compute the PageRank vector and the largest eigenvalue of a left stochastic matrix using the power iteration method, applied to both real and randomly generated web graphs. The implementation uses MPI (Message Passing Interface) for parallelization to explore performance scaling.

---

## 📁 Files in the Repository

- `MPI_GooglePageRank.cpp`: Main C++ file implementing the PageRank algorithm using MPI.
- `Test_GooglePageRank.cpp`: Script for validating correctness using a 6×6 test matrix.
- `script.pbs`: PBS job submission script for running the MPI code on a cluster.
- `README.md`: This file.

---

## 🔧 Compilation Instructions

Compile the main code using the MPI C++ compiler:

```bash
mpic++ -o hpc_code.out MPI_GooglePageRank.cpp
```

---

## 🚀 Running the Code

You can execute the compiled binary as follows:

```bash
mpirun -np <num_procs> ./hpc_code.out <MATRIX_SIZE>
```

Where:
- `<num_procs>` is the number of MPI processes
- `<MATRIX_SIZE>` is the size of the web graph matrix (e.g., 10000)

---

## 📄 PBS Job Script

To run the code on an HPC cluster using PBS:

```bash
qsub script.pbs
```

**PBS Script Overview:**
- Queue: `teachingq`
- Resources: 1 node, 60 CPUs, 248 GB RAM
- Matrix size: 10000×10000
- Wall time: 5 minutes
- Output: Written to `pbs_out.txt`

---

## ✅ Sample Output

**Test Case (6×6 Matrix):**
- Final PageRank vector: `[0.0238, 0.0238, 0.2778, 0.0952, 0.2143, 0.3651]`
- Final Rayleigh Quotient: `1.000000`

**Real Case (10000×10000 Matrix, 60 processors):**
- Final Rayleigh Quotient: `1.000000`
- Execution Time: `0.460119 seconds`

---

## 📈 Performance Scaling

### Strong Scaling (Fixed Problem Size)

| Processors | Execution Time (s) |
|------------|--------------------|
| 1          | 18.04              |
| 2          | 9.03               |
| 4          | 4.75               |
| 12         | 1.91               |
| 60         | 0.46               |

### Weak Scaling (Matrix size ∝ Processors)

| Matrix Size | Processors | Execution Time (s) |
|-------------|------------|--------------------|
| 1291×1291   | 1          | 0.30               |
| 2581×2581   | 4          | 0.33               |
| 10000×10000 | 60         | 0.47               |

---

## 📚 Theoretical Background

The PageRank algorithm is based on the eigenvalue problem of a left stochastic matrix `P` derived from a web link matrix `L`. The power iteration method is used to compute the dominant eigenvector (PageRank vector), satisfying:

```
Pr = λ_max * r
```

Where:
- `r`: PageRank vector
- `λ_max`: Approximate dominant eigenvalue (≈1)

---

## 🧪 Testing

- Random reproducible test cases for `n = 100`, `1000`, `10000` using seeded generation.
- Validation against summation conditions (∑ri = 1, ri ≥ 0).
- Visual inspection of the final `P` matrix and `r` vector for test cases.

---

## 👨‍💻 Author

**Vishal Soni**  
Matriculation Number: 66158  
TU Bergakademie Freiberg

---

## 📜 License

This project is provided as coursework and intended for academic use only.