%% Generalization of the Remez Algorithm
% Approximation space: S = Lin{1, x^2, x^4}
% Find best uniform approximant f_tilde in S for
%   f(x) = exp(3x/5) * cos(2*pi*x/5)
% on [1, 3], tol = 1e-3, kmax = 10.

% --- Test data ---
f    = @(x) exp(3*x/5) .* cos(2*pi*x/5);
a    = 1;
b    = 3;
tol  = 1e-3;
kmax = 10;

% --- Run the generalized Remez algorithm ---
sol = RemezSpecial(f, a, b, tol, kmax);

% --- Display results ---
fprintf('\nBest approximant: f_tilde(x) = c1 + c2*x^2 + c3*x^4\n');
fprintf('  c1 = %.16g\n', sol.coeffs(1));
fprintf('  c2 = %.16g\n', sol.coeffs(2));
fprintf('  c3 = %.16g\n', sol.coeffs(3));
fprintf('Max error (leveling) h = %.16g\n', sol.h);
fprintf('Max error on grid      = %.16g\n', sol.error);
fprintf('Converged in %d iterations\n', sol.iters);

% --- Plot ---
xfine = linspace(a, b, 1000)';
ffine = f(xfine);
pfine = sol.coeffs(1) + sol.coeffs(2)*xfine.^2 + sol.coeffs(3)*xfine.^4;

figure;
subplot(2,1,1); hold on;
plot(xfine, ffine, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x)');
plot(xfine, pfine, 'r--', 'LineWidth', 1.5, 'DisplayName', '\tilde{f}(x)');
legend('show', 'Location', 'best');
title('Function and best uniform approximant');
xlabel('x'); ylabel('y'); hold off;

subplot(2,1,2);
plot(xfine, ffine - pfine, 'b-', 'LineWidth', 1.5);
yline(sol.h, 'r--'); yline(-sol.h, 'r--');
title('Error f(x) - \tilde{f}(x)');
xlabel('x'); ylabel('error');


%% =====================================================================
function sol = RemezSpecial(fun, a, b, tol, kmax)
% REMEZSPECIAL  Generalized Remez algorithm for S = Lin{1, x^2, x^4}.
%   Uses Remez.m / RemezStep.m style approach from tutorials.
%
%   Inputs:
%       fun  – function handle to approximate
%       a, b – interval [a,b] with a > 0
%       tol  – tolerance for convergence
%       kmax – max number of iterations
%
%   Output (struct):
%       sol.coeffs – [c1; c2; c3]  so that f_tilde = c1 + c2*x^2 + c3*x^4
%       sol.h      – final equioscillation error
%       sol.error  – max |f - f_tilde| on fine grid
%       sol.iters  – number of iterations performed

    n    = 3;               % dimension of S: {1, x^2, x^4}
    nref = n + 1;           % need n+1 = 4 reference points

    % --- Fine evaluation grid ---
    M     = 10000;
    xfine = linspace(a, b, M)';
    ffine = fun(xfine);

    % --- Initial reference: Chebyshev-like nodes on [a,b] ---
    j    = 0:(nref - 1);
    refs = (a + b)/2 - (b - a)/2 * cos(j * pi / (nref - 1));
    refs = refs(:);

    coeffs = zeros(n, 1);
    h      = 0;
    iters  = kmax;

    for k = 1:kmax
        % ---- RemezStep: solve leveling equations ----
        [coeffs, h] = RemezStep(fun, refs, n);

        % ---- Evaluate error on fine grid ----
        pfine  = coeffs(1) + coeffs(2)*xfine.^2 + coeffs(3)*xfine.^4;
        efine  = ffine - pfine;
        maxerr = max(abs(efine));

        % ---- Check convergence ----
        if abs(maxerr - abs(h)) < tol
            iters = k;
            break;
        end

        % ---- Exchange: find new reference set ----
        refs = exchange(xfine, efine, nref);
    end

    sol.coeffs = coeffs;
    sol.h      = abs(h);
    sol.error  = maxerr;
    sol.iters  = iters;
end


%% =====================================================================
function [coeffs, h] = RemezStep(fun, refs, n)
% REMEZSTEP  One step of the generalized Remez algorithm.
%   Solves the leveling equations:
%     c1 + c2*x_i^2 + c3*x_i^4 + (-1)^i * h = f(x_i),  i = 1..n+1
%
%   Inputs:
%       fun  – function handle
%       refs – current reference points (n+1 vector)
%       n    – basis dimension (3 for {1, x^2, x^4})
%
%   Outputs:
%       coeffs – [c1; c2; c3]
%       h      – equioscillation error

    nref = length(refs);
    A    = zeros(nref, nref);
    rhs  = fun(refs(:));

    for i = 1:nref
        A(i, 1) = 1;               % basis function u1 = 1
        A(i, 2) = refs(i)^2;       % basis function u2 = x^2
        A(i, 3) = refs(i)^4;       % basis function u3 = x^4
        A(i, 4) = (-1)^i;          % alternating sign for h
    end

    params = A \ rhs;
    coeffs = params(1:n);
    h      = params(n + 1);
end


%% =====================================================================
function refs = exchange(x, e, nref)
% EXCHANGE  Find nref new reference points with alternating extrema.
%   Implements the Remez exchange step.

    N = length(x);

    % --- Collect local extrema indices (including endpoints) ---
    idx = [1];                              %#ok<NBRAK>
    for i = 2:N-1
        if (e(i) > e(i-1) && e(i) > e(i+1)) || ...
           (e(i) < e(i-1) && e(i) < e(i+1))
            idx(end+1, 1) = i;             %#ok<AGROW>
        end
    end
    idx(end+1, 1) = N;
    vals = e(idx);

    % --- Merge consecutive same-sign extrema (keep larger |e|) ---
    fi = idx(1);
    fv = vals(1);
    for i = 2:length(idx)
        if sign(vals(i)) == sign(fv(end))
            if abs(vals(i)) > abs(fv(end))
                fi(end) = idx(i);
                fv(end) = vals(i);
            end
        else
            fi(end+1, 1) = idx(i);         %#ok<AGROW>
            fv(end+1, 1) = vals(i);         %#ok<AGROW>
        end
    end

    % --- Pick nref consecutive alternating extrema maximising min |e| ---
    if length(fi) >= nref
        best_min = -inf;
        best_s   = 1;
        for s = 1:(length(fi) - nref + 1)
            cur_min = min(abs(fv(s:s + nref - 1)));
            if cur_min > best_min
                best_min = cur_min;
                best_s   = s;
            end
        end
        refs = x(fi(best_s:best_s + nref - 1));
    else
        % Fallback: equispaced
        refs = linspace(x(1), x(end), nref)';
    end
end
