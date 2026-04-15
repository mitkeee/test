%% Question 1: Generalized Remez Algorithm
% System of functions S = {1, x^2, x^4} on [a, b] = [1.3, 36/10]
% f(x) = ln(23x/10) * cos(pi*x/2)

% --- Test parameters ---
f   = @(x) log(23*x/10) .* cos(pi*x/2);
a   = 1.3;
b   = 36/10;       % = 3.6
tol = 1e-6;
kmax = 15;

% --- Run the Remez algorithm ---
sol = RemezSpecial(f, a, b, tol, kmax);

% (a) Display coefficients and the constant term
fprintf('\n(a) Constant term c1 = %.17g\n', sol.coeffs(1));

% (b) Maximum absolute approximation error
fprintf('(b) Max absolute error = %.6f\n', sol.error);

% (c) The figure is plotted inside RemezSpecial
fprintf('(c) Figure plotted with f(x) and all iteration polynomials.\n');


%% =====================================================================
function sol = RemezSpecial(fun, a, b, tol, kmax)
% REMEZSPECIAL  Best uniform approximation via the generalized Remez
%               algorithm with Chebyshev system S = {1, x^2, x^4}.
%
%   sol = RemezSpecial(fun, a, b, tol, kmax)
%
%   Inputs:
%       fun  – function handle to approximate
%       a, b – interval endpoints
%       tol  – convergence tolerance (on equioscillation)
%       kmax – maximum number of Remez iterations
%
%   Output (struct):
%       sol.coeffs     – final coefficients [c1; c2; c3]
%       sol.error      – max |f(x) - p(x)| on the fine grid
%       sol.all_coeffs – coefficients at every iteration (for plotting)

    n    = 3;           % basis size: u1=1, u2=x^2, u3=x^4
    nref = n + 1;       % equioscillation requires n+1 = 4 reference points

    % --- Initial reference: Chebyshev nodes mapped to [a, b] ---
    j    = 0:nref-1;
    refs = (a+b)/2 - (b-a)/2 * cos(j*pi/(nref-1));
    refs = refs(:);

    % --- Fine evaluation grid ---
    M     = 10000;
    xfine = linspace(a, b, M)';
    ffine = fun(xfine);

    % --- Storage for every iteration ---
    all_coeffs = zeros(kmax, n);
    num_iters  = kmax;              % will be updated if we converge early

    for k = 1:kmax
        % ---- Step 1: Solve the leveling equations ----
        % c1·1 + c2·xi^2 + c3·xi^4 + (-1)^i · h = f(xi),  i = 1..4
        A   = zeros(nref, nref);
        rhs = fun(refs);
        for i = 1:nref
            A(i,1) = 1;
            A(i,2) = refs(i)^2;
            A(i,3) = refs(i)^4;
            A(i,4) = (-1)^i;       % alternating sign for the error h
        end
        params = A \ rhs;
        coeffs = params(1:n);
        h      = params(n+1);

        all_coeffs(k,:) = coeffs';

        % ---- Step 2: Evaluate error on the fine grid ----
        pfine  = coeffs(1) + coeffs(2)*xfine.^2 + coeffs(3)*xfine.^4;
        efine  = ffine - pfine;
        maxerr = max(abs(efine));

        % ---- Step 3: Check convergence ----
        if abs(maxerr - abs(h)) < tol
            num_iters = k;
            break;
        end

        % ---- Step 4: Exchange – new reference points ----
        refs = newReference(xfine, efine, nref);
    end

    % ---- Plotting: f and the polynomial at EACH iteration ----
    figure; hold on;
    plot(xfine, ffine, 'k-', 'LineWidth', 2, 'DisplayName', 'f(x)');
    cmap = lines(num_iters);        % distinct colour per iteration
    for k = 1:num_iters
        pk = all_coeffs(k,1) + all_coeffs(k,2)*xfine.^2 ...
                              + all_coeffs(k,3)*xfine.^4;
        plot(xfine, pk, 'Color', cmap(k,:), ...
             'DisplayName', sprintf('Iter %d', k));
    end
    legend('show', 'Location', 'best');
    title('Generalized Remez: f(x) and approximating polynomials');
    xlabel('x'); ylabel('y');
    hold off;

    % ---- Return results ----
    sol.coeffs     = coeffs;
    sol.error      = maxerr;
    sol.h          = h;
    sol.all_coeffs = all_coeffs(1:num_iters,:);

    fprintf('\nCoefficients of p(x) = c1 + c2·x^2 + c3·x^4:\n');
    fprintf('  c1 (constant) = %.17g\n', coeffs(1));
    fprintf('  c2 (x^2)      = %.17g\n', coeffs(2));
    fprintf('  c3 (x^4)      = %.17g\n', coeffs(3));
    fprintf('Max |error|     = %.6f\n', maxerr);
    fprintf('Converged at iteration %d / %d\n', num_iters, kmax);
end


%% =====================================================================
function refs = newReference(x, e, nref)
% NEWREFERENCE  Find nref new reference points where |e(x)| is largest
%               with alternating signs (equioscillation exchange step).

    N = length(x);

    % --- Collect ALL local-extremum indices (including endpoints) ---
    idx = 1;                                % left endpoint
    for i = 2:N-1
        if (e(i) > e(i-1) && e(i) > e(i+1)) || ...
           (e(i) < e(i-1) && e(i) < e(i+1))
            idx(end+1,1) = i;              %#ok<AGROW>
        end
    end
    idx(end+1,1) = N;                       % right endpoint
    vals = e(idx);

    % --- Merge consecutive same-sign extrema (keep larger |e|) ---
    fi = idx(1);
    fv = vals(1);
    for i = 2:length(idx)
        if sign(vals(i)) == sign(fv(end))
            % same sign → keep the one with bigger |e|
            if abs(vals(i)) > abs(fv(end))
                fi(end) = idx(i);
                fv(end) = vals(i);
            end
        else
            fi(end+1,1) = idx(i);           %#ok<AGROW>
            fv(end+1,1) = vals(i);           %#ok<AGROW>
        end
    end

    % --- Sliding window: pick nref consecutive alternating extrema ---
    %     that maximise the smallest |e| in the window
    if length(fi) >= nref
        best_min = -inf;
        best_s   = 1;
        for s = 1:length(fi) - nref + 1
            cur_min = min(abs(fv(s:s+nref-1)));
            if cur_min > best_min
                best_min = cur_min;
                best_s   = s;
            end
        end
        refs = x(fi(best_s:best_s+nref-1));
    else
        % fallback: equispaced
        refs = linspace(x(1), x(end), nref)';
    end
end
