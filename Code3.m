%% Question 3: Piecewise Linear Interpolation with Hat Functions
%
% Operator I_n : C([a,b]) -> S_{1,x}
%   I_n f(x) = sum_{i=0}^{n} f(x_i) H_i(x)
% where H_i are hat (tent) basis functions.
% This is piecewise linear interpolation at nodes x_n.
%
% Nodes: x_n = {(2i - n)/n : i = 0,1,...,n}  (equispaced on [-1,1])
% k = 14
% f(x) = |x| * cos(x^2)   on [-1, 1]
% Plot f and I_j f for j = 3, 4, ..., k

% --- Parameters ---
f = @(x) abs(x) .* cos(x.^2);
a = -1;  b = 1;
k = 14;

% --- Fine grid for plotting f ---
xfine = linspace(a, b, 1000)';
ffine = f(xfine);

% --- Plot ---
figure; hold on;
plot(xfine, ffine, 'k-', 'LineWidth', 2.5, 'DisplayName', 'f(x)');

cmap = lines(k - 3 + 1);   % distinct colours for j = 3..k

for j = 3:k
    % Compute the nodes x_n for n = j
    % x_i = (2i - n)/n,  i = 0,1,...,n
    i_vec = 0:j;
    xnodes = (2*i_vec - j) / j;    % equispaced from -1 to 1

    % Evaluate f at nodes
    fnodes = f(xnodes);

    % I_j f is piecewise linear interpolation through (xnodes, fnodes)
    % interp1 with 'linear' does exactly this
    Ijf = interp1(xnodes, fnodes, xfine, 'linear');

    plot(xfine, Ijf, 'Color', cmap(j - 2, :), ...
         'DisplayName', sprintf('I_{%d}f', j));
end

legend('show', 'Location', 'best');
title('f(x) = |x| cos(x^2) and piecewise linear approximations I_j f');
xlabel('x'); ylabel('y');
hold off;

fprintf('Figure plotted: f(x) and I_j f for j = 3, 4, ..., %d\n', k);
fprintf('Enter 1 to attach the figure to the zip file.\n');
