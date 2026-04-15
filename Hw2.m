%% Least Squares Method in the Space of Polynomials
% Inner product: <f,g> = sum_{i=1}^{N} f(x_i)*g(x_i)
% Find best least squares polynomial approximant of degree <= n.

% --- Test data ---
f = @(x) exp(3*x/5) .* cos(2*pi*x/5);
a = 1; b = 4;
x = a:1/5:b;  % x = (1, 6/5, 7/5, ..., 19/5, 4)

% --- Degree 4 ---
sol4 = leastSquares(f, x, 4);
fprintf('Degree 4 coefficients (a0 + a1*x + ... + a4*x^4):\n');
for k = 1:length(sol4)
    fprintf('  a%d = %.16g\n', k-1, sol4(k));
end

% --- Degree 5 ---
sol5 = leastSquares(f, x, 5);
fprintf('\nDegree 5 coefficients (a0 + a1*x + ... + a5*x^5):\n');
for k = 1:length(sol5)
    fprintf('  a%d = %.16g\n', k-1, sol5(k));
end

% --- Plot ---
xfine = linspace(a, b, 1000);
ffine = f(xfine);
p4fine = polyval(flip(sol4), xfine);
p5fine = polyval(flip(sol5), xfine);

figure;
subplot(2,1,1); hold on;
plot(xfine, ffine, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x)');
plot(xfine, p4fine, 'r--', 'LineWidth', 1.5, 'DisplayName', 'degree 4');
plot(xfine, p5fine, 'b-.', 'LineWidth', 1.5, 'DisplayName', 'degree 5');
plot(x, f(x), 'ko', 'MarkerSize', 6, 'DisplayName', 'data points');
legend('show', 'Location', 'best');
title('Function and least squares approximants');
xlabel('x'); ylabel('y'); hold off;

subplot(2,1,2); hold on;
plot(xfine, ffine - p4fine, 'r-', 'LineWidth', 1.5, 'DisplayName', 'error deg 4');
plot(xfine, ffine - p5fine, 'b-', 'LineWidth', 1.5, 'DisplayName', 'error deg 5');
legend('show', 'Location', 'best');
title('Approximation error');
xlabel('x'); ylabel('error'); hold off;


%% =====================================================================
function sol = leastSquares(f, x, n)
% LEASTSQUARES  Best least squares polynomial approximant.
%   sol = leastSquares(f, x, n) returns the (n+1)-vector of coefficients
%   [a0; a1; ...; an] of the polynomial p(x) = a0 + a1*x + ... + an*x^n
%   that minimises ||f - p||^2 = sum_i (f(x_i) - p(x_i))^2.
%
%   Inputs:
%       f – function handle
%       x – vector of sample points
%       n – polynomial degree
%
%   Output:
%       sol – coefficient vector [a0; a1; ...; an]

    x = x(:);           % column vector
    N = length(x);
    fvals = f(x);        % function values at sample points

    % Build Vandermonde matrix V: V(i,j) = x_i^(j-1), j = 1..n+1
    V = zeros(N, n+1);
    for j = 0:n
        V(:, j+1) = x.^j;
    end

    % Solve normal equations: V'*V * sol = V' * fvals
    G   = V' * V;        % Gram matrix
    rhs = V' * fvals;    % right-hand side
    sol = G \ rhs;
end
