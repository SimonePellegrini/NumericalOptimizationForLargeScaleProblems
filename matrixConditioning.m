

clearvars; close all; clc;
rng(361840);

dims = [2, 1000, 10000, 100000]; 



% Number of random samples per dimension
num_samples = 6; 

% Definition of the hessian function
hessian_func = @problem76_hess; 


% Array to store the final average results
avg_cond_results = zeros(length(dims), 1);

for i = 1:length(dims)

    n = dims(i);
    % Center of the hypercube
    center_val = 2*ones(n,1);
    
    % Pre-allocate vector for samples
    current_cond_values = zeros(num_samples, 1);
    
    fprintf('Processing n = %-7d ... ', n);
    t_start = tic; % Start timer for this batch
    
    % Generate the random points
    for k = 1:num_samples
        
        x_sample = center_val -1 + (2 * rand(n, 1));
        
        %  Compute the Hessian
        H = hessian_func(x_sample);
        
        % Compute the conditioning number.
        c = condest(H); 
        current_cond_values(k) = c;
    end
    
    % Calculate average
    avg_cond_results(i) = mean(current_cond_values);
    
    t_end = toc(t_start);
    fprintf('Done (%.2f s). Avg Cond: %.2e\n', t_end, avg_cond_results(i));
end
figure('Color', 'w', 'Name', 'Condition Number vs Dimension', 'Position', [100, 100, 900, 700]);

% Use log-log plot because dimensions and condition numbers
% usually span orders of magnitude.
loglog(dims, avg_cond_results, '-o', ...
    'LineWidth', 3, ...
    'MarkerSize', 12, ...
    'MarkerFaceColor', [0.4660 0.6740 0.1880], ...
    'Color', [0.4660 0.6740 0.1880]);

grid on;
set(gca, 'FontSize', 16, 'LineWidth', 1.5, 'GridAlpha', 0.3, ...
    'GridLineStyle', '-', 'MinorGridLineStyle', 'none');

xlabel('Input Dimension (n)', 'FontSize', 20, 'FontWeight', 'bold');
ylabel('Average Condition Number', 'FontSize', 20, 'FontWeight', 'bold');
title({['Condition Number Scaling']}, ...
    'FontSize', 22, 'FontWeight', 'bold');

y_range = log10(max(avg_cond_results)) - log10(min(avg_cond_results));
offset_factor = 0.06; 

for i = 1:length(dims)
    label_y = avg_cond_results(i) / (10^(y_range * offset_factor));
    
    text(dims(i), label_y, sprintf('%.2e', avg_cond_results(i)), ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'center', ...
        'FontSize', 14, 'FontWeight', 'bold');
end

xlim([min(dims)*0.5, max(dims)*1.8]);
ylim([min(avg_cond_results)*0.2, max(avg_cond_results)*3]);

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [9, 7]);

fprintf('\nAnalysis Complete.\n');