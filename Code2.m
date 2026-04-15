%% Question 2: Discrete Weighted Least Squares Polynomial Approximation
%
% Inner product: <f,g> = sum_{i=1}^{N} rho(x_i) f(x_i) g(x_i)
% Find polynomial p of degree <= n that best approximates f
% in the least-squares sense with respect to this inner product.

% --- Define function and weight ---
f   = @(x) 24 * sin(pi * x / 3).^2;
rho = @(x) 3 * x.^2 / pi;

% --- Parameters ---
a = -1;  b = 1;        % interval [a, b]
n = 4;                  % polynomial degree

% --- Define equispaced data points ---
% x_0 = -1, x_{i+1} = x_i + h, vector x = (x_i)_{i=1}^N
N = 20;
h = (b - a) / N;                  % step size = 0.1
x = (a + h : h : b)';             % x_1, ..., x_N  (excludes x_0 = a)

% --- Compute least squares approximation ---
sol = leastSquaresRho(f, rho, x, n);

% Flip coefficients for polyval: [cn, ..., c1, c0]
pcoeffs = flip(sol.coeffs(:)');

% --- Display all coefficients ---
fprintf('\nCoefficients (c0, c1, ..., c%d):\n', n);
for j = 0:n
    fprintf('  c%d = %.16g\n', j, sol.coeffs(j+1));
end

% (a) Leading coefficient (coefficient of x^n)
fprintf('\n(a) Leading coefficient = %.15g\n', sol.coeffs(n+1));
% Expected: -7.441139663894256

% (b) L2 norm of the error with the discrete weighted inner product
%     ||f - p||_rho = sqrt( sum rho(x_i) * (f(x_i) - p(x_i))^2 )
p_at_x   = polyval(pcoeffs, x);
err_at_x = f(x) - p_at_x;
L2_error = sqrt(sum(rho(x) .* err_at_x.^2));
fprintf('(b) L2 error norm = %.16g\n', L2_error);
% Expected: 0.0670433621606801

% (c) Infinity-norm error ||p - f||_{inf,[-1,1]}
%     approximated numerically on 101 equispaced points
xeval    = a + (0:100)' * (b - a) / 100;   % 101 points on [a, b]
inf_error = max(abs(polyval(pcoeffs, xeval) - f(xeval)));
fprintf('(c) ||p - f||_inf  = %.16g\n', inf_error);
% Expected: 0.1164831513602008


%% =====================================================================
function sol = leastSquaresRho(fun, rho, x, n)
% LEASTSQUARESRHO  Best least-squares polynomial approximation with
%                  respect to a weighted discrete inner product.
%
%   sol = leastSquaresRho(fun, rho, x, n)
%
%   Inputs:
%       fun – function handle to approximate
%       rho – weight function handle  (rho(x_i) >= 0)
%       x   – vector of data points
%       n   – maximum polynomial degree
%
%   Output (struct):
%       sol.coeffs – coefficients [c0; c1; ...; cn] in ascending powers
%                    p(x) = c0 + c1*x + c2*x^2 + ... + cn*x^n

    x = x(:);                       % ensure column vector
    M = length(x);

    % Function values and weights at the nodes
    fvals = fun(x);
    w     = rho(x);

    % Build Vandermonde matrix  V(i,j) = x_i^(j-1),  j = 1..n+1
    V = zeros(M, n + 1);
    for j = 0:n
        V(:, j+1) = x.^j;
    end

    % Weighted least squares:  minimise sum_i w_i (f_i - p_i)^2
    % Multiply rows by sqrt(w) and solve the overdetermined system
    sqrtw  = sqrt(w);
    A      = sqrtw .* V;             % weighted Vandermonde
    b      = sqrtw .* fvals;         % weighted RHS
    coeffs = A \ b;                  % solve via QR factorisation

    sol.coeffs = coeffs;
end
