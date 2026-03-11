function W = compute_W_o(A, mu)

% ----------------------------------------------------------
% Computes matrix W where each column solves
%
%   minimize   (1/n)||w_j||_2^2
%   subject to || (1/n)A'*w_j - e_j ||_inf <= mu
%
% Inputs
%   A  : n x p measurement matrix
%   mu : constraint tolerance
%
% Output
%   W  : n x p matrix
% ----------------------------------------------------------

[n, p] = size(A);

W = zeros(n,p);

for j = 1:p

    e_j = zeros(p,1);
    e_j(j) = 1;

    cvx_begin quiet
        variable w(n)

        minimize( (1/n) * (w' * w) )

        subject to
            norm((1/n) * A' * w - e_j, Inf) <= mu
    cvx_end

    W(:,j) = w;

end

end