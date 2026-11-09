function ExactGradientAndHessian(file_name, settings, funcs)
% Runs experiments with Exact Gradient and Hessian.

rng(361840)

dims = settings.problem_dims;
p_cond_list = settings.use_precond;
n_runs = settings.n_runs;
sp = settings.solver_params;

output_dir = fullfile('.', 'Experiments');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end
f_csv = fullfile(output_dir, file_name + ".csv");
favg_csv = fullfile(output_dir, file_name + "_avg.csv");

% Generate points
points = cell(length(dims), n_runs+1);
for j = 1:length(dims)
    points{j,1} = funcs.x0_gen(dims(j), 0);
    for i = 1:n_runs
        points{j,i+1} = funcs.x0_gen(dims(j), i);
    end
end

f = funcs.obj; g = funcs.grad; h_func = funcs.hess;
methods = settings.methods;

fid = fopen(f_csv, 'w');
fidAvg = fopen(favg_csv,'w');
if fid == -1, error('File open error.'); end

h_raw = 'Method,N,Precond,Time_s,GradientNorm,Iter,ConvRate_Rho,Success';
h_avg = 'Method,N,Precond,AvgTime_s,AvgGradientNorm,AvgIter,AvgConvRate_Rho,SuccessRate';
fprintf(fid, '%s\n', h_raw); fprintf(fidAvg, '%s\n', h_avg);

% Main Loops
for m_idx = 1:length(methods)
    for p_c = p_cond_list
        for j = 1:length(dims)
            n = dims(j);
            m_name = methods{m_idx};
            
            stats = zeros(n_runs, 4); 
            n_succ = 0;
            
            for r = 1:n_runs+1
                if strcmp(m_name, 'modifiedNewton') && ~p_c
                    tic;
                    [gnorm, ~, errS, ~, k, ~] = modifiedNewton(points{j,r}, f, g, h_func, ...
                        sp.max_iter, sp.tolgrad, sp.c1, sp.rho, sp.btmax, sp.beta);
                    t = toc;
                elseif strcmp(m_name, 'truncatedNewton')
                    tic;
                    [gnorm, ~, errS, ~, k, ~] = truncatedNewton(points{j,r}, f, g, h_func, ...
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
                
                fprintf(fid, '%s,%d,%d,%.4f,%.2e,%d,%.4f,%d\n', ...
                    m_name, n, p_c, t, gnorm, k, rho, success);
                fprintf('%d,%d\n',n,k)
            end
            
            if strcmp(m_name, "truncatedNewton") || ~p_c
                if n_succ > 0, means = mean(stats(1:n_succ, :), 1); else, means = zeros(1,4); end
                fprintf(fidAvg, '%s,%d,%d,%.4f,%.2e,%.2f,%.4f,%.2f\n', ...
                    m_name, n, p_c, means(2), means(4), means(1), means(3), n_succ/(n_runs+1));
            end
        end
    end
end
fclose(fid); fclose(fidAvg);
end