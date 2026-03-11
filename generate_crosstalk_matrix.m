function A_tilde = generate_crosstalk_matrix(n, p, zeta)
% generate_crosstalk_matrix: Produces sensing matrix with cross-talk
%
% A is Rademacher (±1 entries)
% Cross-talk model: 
%   a_tilde(j,:) = sqrt(1 - zeta)*a(j,:) + sqrt(zeta)*ones(1,p)
%
% INPUT:
%   n     - number of rows (measurements)
%   p     - number of columns (items/variables)
%   zeta  - cross-talk parameter in [0,1]
%
% OUTPUT:
%   A_tilde - n × p modified sensing matrix

    % Generate Rademacher matrix A (±1)
    A = sign(randn(n,p));   % each entry = +1 or -1

    % Row of ones (cross-talk leakage vector)
    ones_row = ones(1,p);

    % Preallocate
    A_tilde = zeros(n,p);

    % Construct cross-talk contaminated matrix row by row
    for j = 1:n
        A_tilde(j,:) = sqrt(1 - zeta) * A(j,:) + sqrt(zeta) * ones_row;
    end

end