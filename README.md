# Numerical Optimization for Large-Scale Problems

This repository contains MATLAB implementations and numerical experiments developed as part of the *Numerical Optimization for Large-Scale Problems* project.

The project investigates second-order optimization methods, focusing on their numerical performance, convergence behavior, and scalability when dealing with large-dimensional optimization problems.

## Implemented Methods

Two Newton-type optimization methods are considered:

- **Modified Newton Method**: a regularized Newton method that ensures positive definiteness of the Hessian by adding a suitable multiple of the identity matrix. The resulting linear system is solved using Cholesky factorization.
- **Truncated Newton Method**: an iterative Newton-type method based on an approximate solution of the Newton system, using the Preconditioned Conjugate Gradient (PCG) method.

Both methods employ an Armijo backtracking line search to ensure sufficient decrease of the objective function.

## Gradient and Hessian Approximations

The project investigates the impact of approximating first- and second-order derivatives on optimization performance.

Three experimental configurations are considered:

1. **Exact Gradient and Exact Hessian**: both derivatives are computed analytically.
2. **Approximate Gradient and Approximate Hessian**: derivatives are estimated using finite differences.
3. **Exact Gradient and Approximate Hessian**: the gradient is computed analytically, while the Hessian is approximated numerically.

Different finite-difference step sizes and relative/absolute approximation strategies are investigated.

## Preconditioning

The Truncated Newton method is also evaluated with and without preconditioning.

The experiments aim to investigate how preconditioning affects the convergence behavior and computational efficiency of iterative linear-system solvers, particularly as the problem dimension increases.

## Numerical Experiments

The repository includes experimental routines to evaluate and compare the implemented methods.

The main performance indicators include:

- **Execution Time**: computational time required by each method.
- **Gradient Norm**: norm of the gradient at the final iterate.
- **Number of Iterations**: iterations required to satisfy the stopping criterion.
- **Convergence Rate**: empirical convergence rate estimated from the sequence of optimization errors.
- **Success Rate**: fraction of runs that converge before reaching the maximum iteration limit.

Experiments are performed using multiple initial points and different problem dimensions. Individual results and aggregated statistics are exported to CSV files.

### Visualization

Contour plots of the objective function and corresponding optimization trajectories can be generated to visualize the behavior of the methods from different starting points.

### Matrix Conditioning

An additional experiment investigates how the condition number of the Hessian matrix changes with the problem dimension, providing insight into the numerical properties of the underlying optimization problem.

## Repository Structure

| File | Description |
|---|---|
| `modifiedNewton.m` | Implementation of the Modified Newton method. |
| `ExactGradientAndHessian.m` | Experiments using exact first- and second-order derivatives. |
| `ExactGradientApxHessian.m` | Experiments using an exact gradient and an approximate Hessian. |
| `ApxGradientAndHessian.m` | Experiments using approximate gradients and Hessians. |
| `genGraphs.m` | Generates contour plots and optimization trajectories. |
| `graphsProblem76.m` | Configuration and execution of visualization experiments for Problem 76. |
| `matrixConditioning.m` | Studies Hessian conditioning across different dimensions. |

The experimental routines save their output in an `Experiments/` directory, while visualization results are saved in the configured output folder.

## Requirements

- MATLAB
- Optimization problem definitions and derivative functions.
- Implementations of the Modified Newton and Truncated Newton methods.

## Getting Started

Clone the repository:

```bash
git clone https://github.com/SimonePellegrini/NumericalOptimizationForLargeScaleProblems.git
```

Navigate to the project directory:

```bash
cd NumericalOptimizationForLargeScaleProblems
```

Run the desired experimental driver after configuring the objective function, its derivatives, and the corresponding solver parameters.

## Author

**Simone Pellegrini**

## Purpose

This project explores the practical performance of second-order numerical optimization methods, examining the trade-offs between derivative approximation, computational efficiency, preconditioning, and convergence behavior in large-scale optimization problems.
