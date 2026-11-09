function genGraphs(file_prefix, settings, funcs, plot_opts)
% Generates contour plots with optimization paths.

rng(361840);

% Extract settings
methods_list = settings.methods;
use_precond  = settings.use_precond;
h_val        = settings.h_val;
relative     = settings.relative;
apxOnlyHess  = settings.apxOnlyHess;
n_starts     = settings.n_starts;
sp           = settings.solver_params;

% Setup Output Directory
out_dir = fullfile('.', plot_opts.folder_name);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Generate Grid for Contours
[xg, yg] = meshgrid(linspace(plot_opts.x_lim(1), plot_opts.x_lim(2), 100), ...
                    linspace(plot_opts.y_lim(1), plot_opts.y_lim(2), 100));
zg = arrayfun(@(x,y) funcs.obj([x; y]), xg, yg);

% Generate Starting Points (Fixed for all loops)
X0 = zeros(2, n_starts + 1);
X0(:, 1) = funcs.x0_gen(0); % Center/Reference
for i = 1:n_starts
    X0(:, i+1) = funcs.x0_gen(i);
end

colors = plot_opts.colors;

% Main Loops
for m = 1:length(methods_list)
    for p = 1:length(use_precond)
        for hh = 1:length(h_val)
            for r = 1:length(relative)
                for a = 1:length(apxOnlyHess)
                    
                    methodName = methods_list(m);
                    isPrecond  = use_precond(p);
                    h          = h_val(hh);
                    isRel      = relative(r);
                    isApxH     = apxOnlyHess(a);
                    
                    % Skip invalid combinations
                    if strcmp(methodName, "modifiedNewton") && isPrecond, continue; end
                    
                    % Select Gradient and Hessian functions based on flags
                    f = funcs.obj;
                    if h == 0
                        % Case: Exact
                        gradf    = funcs.grad_exact;
                        hessianf = funcs.hess_exact;
                    else
                        if isApxH
                            % Case: Hybrid (Exact Grad, Apx Hess)
                            gradf    = funcs.grad_exact;
                            hessianf = @(x) funcs.hess_gen(x, h, isRel);
                        else
                            % Case: Fully Approx (FD Grad, Static FD Hess)
                            gradf    = @(x) funcs.grad_apx(x, h, isRel);
                            hessianf = funcs.hess_static; 
                        end
                    end
                    
                    % Plotting setup
                    fig = figure('Visible', 'off');
                    hold on; grid on;
                    set(gca, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
                    axis tight;
                    
                    contour(xg, yg, zg, 40, 'LineColor', [0.65 0.85 0.65], 'LineWidth', 0.8);
                    
                    legend_entries = cell(n_starts + 1, 1);
                    h_paths = gobjects(n_starts + 1, 1);
                    
                    % Run solver for each starting point
                    for i = 1:n_starts + 1
                        x0 = X0(:, i);
                        
                        if strcmp(methodName, "modifiedNewton")
                            [~, ~, ~, ~, k, xseq] = modifiedNewton( ...
                                x0, f, gradf, hessianf, ...
                                sp.max_iter, sp.tolgrad, sp.c1, sp.rho, sp.btmax, sp.beta);
                        else
                            [~, ~, ~, ~, k, xseq] = truncatedNewton( ...
                                x0, f, gradf, hessianf, ...
                                sp.max_iter, sp.tolgrad, sp.c1, sp.rho, sp.btmax, sp.pcg_max_iter, isPrecond);
                        end
                        
                        % Draw Path
                        col = colors(mod(i-1, size(colors,1))+1, :);
                        h_paths(i) = plot(xseq(:, 1), xseq(:, 2), '-o', ...
                            'Color', col, 'MarkerFaceColor', col, ...
                            'LineWidth', 1.4, 'MarkerSize', 4);
                        
                        plot(xseq(1, 1), xseq(1, 2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
                        plot(xseq(end, 1), xseq(end, 2), 'r*', 'MarkerSize', 8);
                        
                        legend_entries{i} = sprintf('Start (%.2f, %.2f), Iter: %d', x0(1), x0(2), k);
                    end
                    
                    lgd = legend(h_paths, legend_entries);
                    set(lgd, 'Location', 'northwest', 'FontSize', 12);
                    
                    % Filename construction
                    fname = sprintf('%s_%s_precond_%d_h_%g_rel_%d_apxH_%d.pdf', ...
                        file_prefix, methodName, isPrecond, h, isRel, isApxH);
                    
                    exportgraphics(fig, fullfile(out_dir, fname), 'ContentType', 'vector');
                    close(fig);
                end
            end
        end
    end
end
fprintf('All contour plots generated in %s\n', out_dir);
end