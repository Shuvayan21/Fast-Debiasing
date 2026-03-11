function W = compute_W_e(A)
% ----------------------------------------------------------
% Computes W using the exact formula:
%
%   w_j = [ n*(1 - mu) / ||a_j||_2^2 ] * a_j
%
% where mu = rho/(1+rho), and
%
%   rho = max_{i != j} |a_i' a_j| / ||a_j||_2^2.
%
% Inputs:
%   A : n x p matrix
%
% Output:
%   W : n x p matrix
% ----------------------------------------------------------

[n, p] = size(A);

% Compute rho(A)
rho = compute_rho(A);

% mu = rho / (1 + rho)
mu = rho / (1 + rho);

% Prefactor = n*(1 - mu)
pref = n * (1 - mu);

% Precompute denominator ||a_j||_2^2
norm_cols_sq = sum(A.^2, 1);

% Allocate W
W = zeros(n, p);

for j = 1:p
    W(:, j) = (pref / norm_cols_sq(j)) * A(:, j);
end

end