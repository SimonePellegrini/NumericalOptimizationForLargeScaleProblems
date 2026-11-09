% Visual Experiment Driver
clear; clc; close all;
rng(361840);

%% 1. Experiment Settings
settings = struct();
settings.methods     = ["modifiedNewton", "truncatedNewton"];
settings.use_precond = [false, true];
settings.h_val       = [0, 1e-4, 1e-8, 1e-12]; 
settings.relative    = [true, false];
settings.apxOnlyHess = [true, false];
settings.n_starts    = 5;

% Solver Parameters
settings.solver_params = struct(...
    'max_iter', 1000, 'tolgrad', 1e-6, 'c1', 0.1, ...
    'rho', 0.5, 'btmax', 50, 'pcg_max_iter', 100, 'beta', 0.5);

%% 2. Plot Options
plot_opts = struct();
plot_opts.folder_name = 'trigonometric_Images';
plot_opts.x_lim       = [-5, 3];
plot_opts.y_lim       = [-3, 3];
plot_opts.colors      = [
    0.93 0.69 0.13; 0.30 0.75 0.93; 0.00 0.60 0.50;
    0.90 0.40 0.70; 0.60 0.60 0.60; 0.20 0.20 0.20
];

%% 3. Function Handles 
funcs = struct();
funcs.obj         = @(x) bandedTrigonometric_func(x);
funcs.grad_exact  = @(x) bandedTrigonometric_grad(x);
funcs.hess_exact  = @(x) bandedTrigonometric_hess(x);

% Approximations handles
funcs.grad_apx    = @(x, h, rel) bandedTrigonometric_FDGrad(x, h, rel);
funcs.hess_static = @(x) bandedTrigonometric_FDHess(x, 1e-4); % Used when fully approx
funcs.hess_gen    = @(x, h, rel) bandedTrigonometric_hessGivenGrad(x, h, rel); % Used when Hybrid

% Start Point Generator
funcs.x0_gen      = @(idx) generate_points_visual(idx);

%% 4. Run
genGraphs("tridiagonal", settings, funcs, plot_opts);

%% Local Helper
function x = generate_points_visual(idx)
    x_center = ones(2, 1);
    if idx == 0
        x = x_center;
    else
        % Logic: x_center + (2 * rand - 1) -> Range [x-1, x+1]
        x = x_center + (2 * rand(2, 1) - 1);
    end
end