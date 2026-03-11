function Sigma = makeSigma(p,xi)
    % p = dimension of the covariance matrix
    Sigma = zeros(p,p);

    for j = 1:p
        for k = 1:p
            if k == j
                Sigma(j,k) = 1;
            elseif any(k == mod(j + (1:5) - 1, p) + 1)  % wrap-around neighbors
                Sigma(j,k) = xi;
            else
                Sigma(j,k) = 0;
            end
        end
    end
end
