% Plot profit-taking ratio decisions of three VPPs over iterations.
clc;

if ~exist('storage', 'var')
    Gauss_Iteration;
end

valid_rows = find(any(storage ~= 0, 2));
storage_plot = storage(valid_rows, :);
iteration_count = valid_rows;

figure('Color', 'w');
plot(iteration_count, storage_plot(:, 1), '-o', ...
    'LineWidth', 2.2, 'MarkerSize', 6);
hold on;
plot(iteration_count, storage_plot(:, 2), '-s', ...
    'LineWidth', 2.2, 'MarkerSize', 6);
plot(iteration_count, storage_plot(:, 3), '-^', ...
    'LineWidth', 2.2, 'MarkerSize', 6);
hold off;

set(gca, 'FontName', 'Calibri', 'FontSize', 12);
xlabel('Iteration count', 'FontName', 'Calibri', 'FontSize', 13);
ylabel('profit-taking ratio', 'FontName', 'Calibri', 'FontSize', 13);
legend({'VPP1', 'VPP2', 'VPP3'}, ...
    'FontName', 'Calibri', 'FontSize', 12, 'Location', 'best');
grid on;
box on;

exportgraphics(gcf, 'profit_taking_ratio_iteration.png', 'Resolution', 300);
savefig(gcf, 'profit_taking_ratio_iteration.fig');
