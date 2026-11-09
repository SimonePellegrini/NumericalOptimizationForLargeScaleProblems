% Master Driver Script
clear; clc; close all;
rng(361840);

% Shared Settings
settings = struct();
settings.problem_dims = [2, 1e3, 1e4, 1e5]; 
settings.n_runs       = 5;
settings.use_precond  = [true, false];
settings.h_val        = [1e-4, 1e-8, 1e-12];
settings.relative     = [false, true]; 
settings.methods      = {"modifiedNewton","truncatedNewton"};

% Function Handles
f_obj      = @(x) problem76_func(x);
grad_exact = @(x) problem76_grad(x);
hess_exact = @(x) problem76_hess(x);
grad_apx   = @(x, h, rel) problem76_FDGrad(x, h, rel);
hess_stat  = @(x) problem76_FDHess(x, 1e-4); 
hess_gen   = @(x, h, rel) problem76_hessianGivenGrad(x, h, rel);

%% 1. Approx Gradient & Hessian
display("0/3");
settings.solver_params = struct('tolgrad', 1e-6, 'max_iter', 1000, 'c1', 1e-4, ...
    'rho', 0.5, 'btmax', 100, 'beta', 1e-3, 'pcg_max_iter', 100);

funcs1 = struct();
funcs1.obj = f_obj; funcs1.grad = grad_apx; funcs1.hess = hess_stat;
funcs1.x0_gen = @(n, k) gen_points_problem76(n, k);

ApxGradientAndHessian("Problem76", settings, funcs1);

%% 2. Exact Gradient & Hessian
display("1/3");
funcs2 = struct();
funcs2.obj = f_obj; funcs2.grad = grad_exact; funcs2.hess = hess_exact;
funcs2.x0_gen = @(n, k) gen_points_problem76(n, k); 

ExactGradientAndHessian("Prolem76", settings, funcs2);

%% 3. Hybrid (Exact Grad / Approx Hess)
display("2/3");
funcs3 = struct();
funcs3.obj = f_obj; funcs3.grad = grad_exact; funcs3.hess_gen = hess_gen;
funcs3.x0_gen = @(n, k) gen_points_problem76(n, k);

ExactGradientApxHessian("Problem76", settings, funcs3);

%% Points generator
display("3/3");
function x = gen_points_problem76(n, run_idx)
    x_ref = 2 * ones(n, 1);
    if run_idx > 0
        x = x_ref -1 + 2 * rand(n, 1);
    else
        x = x_ref;
    end
   
end

