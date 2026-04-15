%% Orthogonal Polynomials
% Gram-Schmidt orthogonalization with discrete inner product
%   <f,g> = sum_{i=1}^{N} f(x_i)*g(x_i)

% --- Test data ---
x = 0:pi/4:2*pi;  % x = (0, pi/4, pi/2, ..., 7*pi/4, 2*pi)
n = 7;             % first 7 orthogonal polynomials

% --- Compute orthogonal polynomials ---
P = orthogonal(x, n);

% --- Display results ---
fprintf('First %d orthogonal polynomials on x = 0:pi/4:2*pi\n\n', n);
for k = 1:n
    fprintf('p_%d(x) = ', k-1);
    coeffs = P(k, :);
    deg = length(coeffs) - 1;
    first = true;
    for j = deg:-1:0
        c = coeffs(deg - j + 1);
        if abs(c) < 1e-12
            continue;
        end
        if first
            if j == 0
                fprintf('%.6g', c);
            elseif j == 1
                fprintf('%.6g*x', c);
            else
                fprintf('%.6g*x^%d', c, j);
            end
            first = false;
        else
            if c > 0
                fprintf(' + ');
            else
                fprintf(' - ');
                c = -c;
            end
            if j == 0
                fprintf('%.6g', c);
            elseif j == 1
                fprintf('%.6g*x', c);
            else
                fprintf('%.6g*x^%d', c, j);
            end
        end
    end
    if first
        fprintf('0');
    end
    fprintf('\n');
end

% --- Verify orthogonality ---
fprintf('\nOrthogonality check (<p_i, p_j> for i ~= j should be ~0):\n');
for i = 1:n
    pi_vals = polyval(P(i,:), x);
    for j = i+1:n
        pj_vals = polyval(P(j,:), x);
        ip = sum(pi_vals .* pj_vals);
        fprintf('  <p_%d, p_%d> = %.2e\n', i-1, j-1, ip);
    end
end


%% =====================================================================
function P = orthogonal(x, n)
% ORTHOGONAL  Gram-Schmidt orthogonal polynomials w.r.t. discrete inner product.
%   P = orthogonal(x, n) returns an n-by-(n) matrix P where the k-th row
%   contains the coefficients of the (k-1)-th orthogonal polynomial in
%   descending power order (suitable for polyval).
%
%   Inner product: <f,g> = sum_{i=1}^{N} f(x_i)*g(x_i)
%
%   Inputs:
%       x – vector of sample points
%       n – number of orthogonal polynomials to compute (p_0, ..., p_{n-1})
%
%   Output:
%       P – n-by-n matrix, row k = coefficients of p_{k-1}(x) in
%           descending power order

    x = x(:)';
    N = length(x);

    % Store polynomials as coefficient vectors (descending powers)
    % and their values at the sample points
    P    = zeros(n, n);      % output matrix
    Vals = zeros(n, N);      % Vals(k,:) = p_{k-1} evaluated at x

    for k = 1:n
        % Start with monomial x^(k-1)
        % In descending power order: [1, 0, 0, ..., 0] with k coefficients
        mono_coeffs = zeros(1, k);
        mono_coeffs(1) = 1;   % x^(k-1)
        mono_vals = x.^(k-1);

        % Subtract projections onto previous orthogonal polynomials
        proj_vals = mono_vals;
        proj_coeffs = zeros(1, n);  % padded to length n (descending)
        % Place mono_coeffs into proj_coeffs (right-aligned, descending)
        proj_coeffs(n-k+1:n) = mono_coeffs;

        for j = 1:k-1
            pj_vals = Vals(j, :);
            % <mono, p_j> / <p_j, p_j>
            coeff = sum(mono_vals .* pj_vals) / sum(pj_vals .* pj_vals);
            proj_vals = proj_vals - coeff * pj_vals;

            % Subtract from coefficient representation
            pj_coeffs = zeros(1, n);
            pj_coeffs(n-j+1:n) = P(j, n-j+1:n);  % already stored
            proj_coeffs = proj_coeffs - coeff * pj_coeffs;
        end

        % Store
        P(k, :) = proj_coeffs;
        Vals(k, :) = proj_vals;
    end
end
