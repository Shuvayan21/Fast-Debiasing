function rho = compute_rho(A)
%COMPUTE_RHO  Computes ρ(A) = max_{i ≠ j} |a_i' * a_j| / ||a_j||_2^2

    % Compute column norms squared
    col_norm_sq = sum(A.^2, 1);  % 1 × p vector of ||a_j||_2^2

    % Compute Gram matrix (inner products between all column pairs)
    G = A' * A;  % p × p matrix, where G(i,j) = a_i' * a_j

    % Normalize each column j by ||a_j||_2^2 in the denominator
    ratio_matrix = abs(G) ./ col_norm_sq;  % divides each column j appropriately

    % Ignore diagonal elements (i = j)
    ratio_matrix(1:size(A,2)+1:end) = 0;

    % Take the maximum over all i ≠ j
    rho = max(ratio_matrix(:));
end
