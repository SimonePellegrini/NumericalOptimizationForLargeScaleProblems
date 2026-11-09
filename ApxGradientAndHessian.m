function ApxGradientAndHessian(file_name, settings, funcs)
% Input:
% - file_name: Base string for the output files
% - settings: Struct containing experiment parameters and solver settings
% - funcs: Struct containing function handles (f, grad, hess, x0)

rng(361840)

% Experiment Parameters
problem_dims = settings.problem_dims;
use_precond  = settings.use_precond;
h_val        = settings.h_val;
relative     = settings.relative;
n_runs       = settings.n_runs;

% Solver Parameters (extracted for easier access)
p = settings.solver_params; 

% Filenames setup
if ~exist('Experiments', 'dir')
    mkdir('Experiments');
end
filename_mod      = fullfile('.', 'Experiments', "apx_mod_" + file_name + ".csv");
filenameavg_mod   = fullfile('.', 'Experiments', "apx_mod_" + file_name + "_avg.csv");
filename_trunc    = fullfile('.', 'Experiments', "apx_trunc_" + file_name + ".csv");
filenameavg_trunc = fullfile('.', 'Experiments', "apx_trunc_" + file_name + "_avg.csv");

% generates the random points used to test the function
points = cell(length(problem_dims), n_runs+1);
for j = 1:length(problem_dims)
    n = problem_dims(j);
    points{j,1} = funcs.x0_gen(n, 0); 
    for i = 1:n_runs
        points{j,i+1} = funcs.x0_gen(n, i);
    end
end

f = funcs.obj; 
apxHessianf = funcs.hess; 
methods = settings.methods;

% File opening
fid_mod = fopen(filename_mod, 'w');
fidAvg_mod = fopen(filenameavg_mod, 'w');
fid_trunc = fopen(filename_trunc, 'w');
fidAvg_trunc = fopen(filenameavg_trunc, 'w');

if fid_mod == -1 || fid_trunc == -1
    error('Impossible to open files for writing.');
end

header_raw = 'Method,N,h,relative,Precond,Time_s,GradientNorm,Iter,ConvRate_Rho,Success';
header_avg = 'Method,N,h,relative,Precond,AvgTime_s,AvgGradientNorm,AvgIter,AvgConvRate_Rho,SuccessRate';
fprintf(fid_mod, '%s\n', header_raw); fprintf(fidAvg_mod, '%s\n', header_avg);
fprintf(fid_trunc, '%s\n', header_raw); fprintf(fidAvg_trunc, '%s\n', header_avg);

for rel_cond = relative
    for p_cond = use_precond
        for j = 1:length(problem_dims)
            for h = h_val
                apxGradf = @(x) funcs.grad(x, h, rel_cond);
                n = problem_dims(j);
                for m_idx = 1:length(methods)
                    method_name = methods{m_idx};
                    if strcmp(method_name, 'modifiedNewton')
                        fid = fid_mod; fidAvg = fidAvg_mod;
                    else
                        fid = fid_trunc; fidAvg = fidAvg_trunc;
                    end
                    
                    avgK = zeros(n_runs,1); avgElapsed = zeros(n_runs,1);
                    avgConv_rate = zeros(n_runs,1); avgGradientNorm = zeros(n_runs,1);
                    n_success = 0;
                    
                    for r = 1:n_runs+1
                        if strcmp(method_name, 'modifiedNewton') && ~p_cond
                            t_start = tic;
                            [gradfk_norm,fk,errorSeq,x_k,k,xseq] = modifiedNewton(points{j,r}, f, apxGradf, apxHessianf, ...
                                p.max_iter, p.tolgrad, p.c1, p.rho, p.btmax, p.beta);
                            elapsed = toc(t_start);
                        elseif strcmp(method_name, 'truncatedNewton')
                            t_start = tic;
                            [gradfk_norm,fk,errorSeq,x_k,k,xseq] = truncatedNewton(points{j,r}, f, apxGradf, apxHessianf, ...
                                p.max_iter, p.tolgrad, p.c1, p.rho, p.btmax, p.pcg_max_iter, p_cond);
                            elapsed = toc(t_start);
                        else
                            continue;
                        end
                        
                        conv_rate_seq = zeros(k-2,1);
                        for i = 2 : k-1
                            if errorSeq(i) > 1e-16 && errorSeq(i-1) > 1e-16
                                conv_rate_seq(i-1) = log(errorSeq(i+1)/errorSeq(i))/log(errorSeq(i)/errorSeq(i-1));
                            end
                        end
                        curr_rho = 0;
                        if ~isempty(conv_rate_seq), curr_rho = conv_rate_seq(end); end
                        
                        if k < p.max_iter
                            success = 1; n_success = n_success + 1;
                            avgConv_rate(n_success) = curr_rho; avgK(n_success) = k;
                            avgElapsed(n_success) = elapsed; avgGradientNorm(n_success) = gradfk_norm;
                        else
                            success = 0;
                        end
                        fprintf('%f,%d,%d\n',h,n,k)
                        fprintf(fid, '%s,%d,%.1e,%d,%d,%.4f,%.2e,%d,%.4f,%d\n', ...
                            method_name, n, h, rel_cond, p_cond, elapsed, gradfk_norm, k, curr_rho, success);
                    end
                    
                    if n_success > 0
                        avgConv_rate = avgConv_rate(1:n_success); avgK = avgK(1:n_success);
                        avgElapsed = avgElapsed(1:n_success);
                    end
                    if strcmp(method_name, "truncatedNewton") || ~p_cond
                      
                        fprintf(fidAvg, '%s,%d,%.1e,%d,%d,%.5f,%.2e,%.2f,%.4f,%.2f\n', ...
                            method_name, n, h, rel_cond, p_cond, mean(avgElapsed), ...
                            mean(avgGradientNorm), mean(avgK), mean(avgConv_rate), n_success/(n_runs+1));
                    end
                end
            end
        end
    end
end
fclose(fid_mod); fclose(fidAvg_mod); fclose(fid_trunc); fclose(fidAvg_trunc);
fprintf('All the experiments are completed!\n');
end