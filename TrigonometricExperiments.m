% Master Driver Script
clear; clc; close all;
rng(361840);

% Shared Settings
settings = struct();
settings.problem_dims = [2, 1e3, 1e4, 1e5]; 
settings.n_runs       = 6;
settings.use_precond  = [true,false];
settings.h_val        = [1e-4, 1e-8, 1e-12];
settings.relative     = [false, true]; 
settings.methods = {'modifiedNewton'};

% Function Handles
f_obj      = @(x) bandedTrigonometric_func(x);
grad_exact = @(x) bandedTrigonometric_grad(x);
hess_exact = @(x) bandedTrigonometric_hess(x);
grad_apx   = @(x, h, rel) bandedTrigonometric_FDGrad(x, h, rel);
hess_apx  = @(x) bandedTrigonometric_FDHess(x, 1e-4); 
hess_gen   = @(x, h, rel) bandedTrigonometric_hessGivenGrad(x, h, rel);
settings.solver_params = struct('tolgrad', 1e-6, 'max_iter', 1000, 'c1', 0.1, ...
    'rho', 0.5, 'btmax',50, 'beta', 0.5, 'pcg_max_iter', 100);

%% 1. Exact Gradient & Hessian

display("0/3")
funcs2 = struct();
funcs2.obj = f_obj; funcs2.grad = grad_exact; funcs2.hess = hess_exact;
funcs2.x0_gen = @(n, k) gen_points_Trigonometric(n, k); 

ExactGradientAndHessian("Trigonometric", settings, funcs2);

%% 2. Approx Gradient & Hessian
display("1/3")
funcs1 = struct();
funcs1.obj = f_obj; funcs1.grad = grad_apx; funcs1.hess = hess_apx;
funcs1.x0_gen = @(n, k) gen_points_Trigonometric(n, k);

ApxGradientAndHessian("Trigonometric", settings, funcs1);


display("2/3")
%% 3. Hybrid (Exact Grad / Approx Hess)

funcs3 = struct();
funcs3.obj = f_obj; funcs3.grad = grad_exact; funcs3.hess_gen = hess_gen;
funcs3.x0_gen = @(n, k) gen_points_Trigonometric(n, k);
display("3/3")
ExactGradientApxHessian("Trigonometric", settings, funcs3);

%% Points generator
function x = gen_points_Trigonometric(n, run_idx)
    x_ref = ones(n, 1);
    if run_idx > 0
        x = x_ref -1 + 2 * rand(n, 1);
    else
        x = x_ref;
    end
end

