function ExactGradientApxHessian(file_name, settings, funcs)
% Runs experiments with Exact Gradient and Approximate Hessian.

rng(361840)

dims = settings.problem_dims;
p_cond_list = settings.use_precond;
h_vals = settings.h_val;
rel_list = settings.relative;
n_runs = settings.n_runs;
sp = settings.solver_params;

output_dir = fullfile('.', 'Experiments');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
f_mod      = fullfile(output_dir, "apxOnlyHessian_modifiedNewton_" + file_name + ".csv");
favg_mod   = fullfile(output_dir, "apxOnlyHessian_modifiedNewton_" + file_name + "avg.csv");
f_trunc    = fullfile(output_dir, "apxOnlyHessian_truncatedNewton_" + file_name + ".csv");
favg_trunc = fullfile(output_dir, "apxOnlyHessian_truncatedNewton_" + file_name + "avg.csv");

% Generate points
points = cell(length(dims), n_runs+1);
for j = 1:length(dims)
    points{j,1} = funcs.x0_gen(dims(j), 0);
    for i = 1:n_runs
        points{j,i+1} = funcs.x0_gen(dims(j), i);
    end
end

f = funcs.obj; g = funcs.grad;
methods = settings.methods;

fid_mod = fopen(f_mod, 'w'); fidAvg_mod = fopen(favg_mod, 'w');
fid_trunc = fopen(f_trunc, 'w'); fidAvg_trunc = fopen(favg_trunc, 'w');
if any([fid_mod, fid_trunc] == -1), error('File open error.'); end

h_raw = 'Method,N,h,relative,Precond,Time_s,gradientNorm,Iter,ConvRate_Rho,Success';
h_avg = 'Method,N,h,relative,Precond,AvgTime_s,avgGradientNorm,AvgIter,AvgConvRate_Rho,SuccessRate';
fprintf(fid_mod, '%s\n', h_raw); fprintf(fidAvg_mod, '%s\n', h_avg);
fprintf(fid_trunc, '%s\n', h_raw); fprintf(fidAvg_trunc, '%s\n', h_avg);

% Main Loops
for m_idx = 1:length(methods)
    for rel = rel_list
        for p_c = p_cond_list
            for j = 1:length(dims)
                for h = h_vals
                    n = dims(j);
                    hess_handle = @(x) funcs.hess_gen(x, h, rel);
                    m_name = methods{m_idx};
   
                    if strcmp(m_name, 'modifiedNewton')
                        fid = fid_mod; fidA = fidAvg_mod;
                    else
                        fid = fid_trunc; fidA = fidAvg_trunc;
                    end
                  
                    stats = zeros(n_runs, 4);
                    n_succ = 0;
                    
                    for r = 1:n_runs+1
                        if strcmp(m_name, 'modifiedNewton') && ~p_c
                            tic;
                            [gnorm, ~, errS, ~, k, ~] = modifiedNewton(points{j,r}, f, g, hess_handle, ...
                                sp.max_iter, sp.tolgrad, sp.c1, sp.rho, sp.btmax, sp.beta);
                            t = toc;
                        elseif strcmp(m_name, 'truncatedNewton')
                            tic;
                            [gnorm, ~, errS, ~, k, ~] = truncatedNewton(points{j,r}, f, g, hess_handle, ...
                                sp.max_iter, sp.tolgrad, sp.c1, sp.rho, sp.btmax, sp.pcg_max_iter, p_c);
                            t = toc;
                        else
                            continue;
                        end
                        
                        rho = 0;
                        if k > 2
                            rhos = log(errS(3:k)./errS(2:k-1)) ./ log(errS(2:k-1)./errS(1:k-2));
                            if ~isempty(rhos), rho = rhos(end); end
                        end
                        
                        success = (k < sp.max_iter);
                        if success
                            n_succ = n_succ + 1;
                            stats(n_succ, :) = [k, t, rho, gnorm];
                        end
                        
                        fprintf(fid, '%s,%d,%.1e,%d,%d,%.4f,%.2e,%d,%.4f,%d\n', ...
                            m_name, n, h, rel, p_c, t, gnorm, k, rho, success);
                    end
                    
                    if strcmp(m_name, "truncatedNewton") || ~p_c
                        if n_succ > 0, means = mean(stats(1:n_succ, :), 1); else, means = zeros(1,4); end
                        fprintf(fidA, '%s,%d,%.1e,%d,%d,%.4f,%.2e,%.2f,%.4f,%.2f\n', ...
                            m_name, n, h, rel, p_c, means(2), means(4), means(1), means(3), n_succ/(n_runs+1));
                    end
                end
            end
        end
    end
end
fclose(fid_mod); fclose(fidAvg_mod); fclose(fid_trunc); fclose(fidAvg_trunc);
end