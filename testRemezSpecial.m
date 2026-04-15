%% Test script for RemezSpecial
% System of functions S = {1, x^2, x^4} on [a, b] = [1.3, 36/10]
% f(x) = ln(23x/10) * cos(pi*x/2)

% --- Test parameters ---
f    = @(x) log(23*x/10) .* cos(pi*x/2);
a    = 1.3;
b    = 36/10;       % = 3.6
tol  = 1e-6;
kmax = 15;

% --- Run the Remez algorithm ---
sol = RemezSpecial(f, a, b, tol, kmax);

% (a) Display coefficients and the constant term
fprintf('\n(a) Constant term c0 = %.14f\n', sol.coeffs(1));

% (b) Maximum absolute approximation error
fprintf('(b) Max absolute error = %.6f\n', sol.error);

% (c) The figure is plotted inside RemezSpecial
fprintf('(c) Figure plotted with f(x) and all iteration polynomials.\n');
fprintf('Enter 1 to attach the figure.\n');
