rng(1)
n_values = 40:10:100; % n = 40, 50, 60, 70, 80, 90, 100
p = 100; % Number of columns
iterations = 1; % Set iterations to 1
 mu_values = linspace(0.2, 0.6, 41); % Various values of mu

% Initialize arrays to store Frobenius norm ratios and thresholds
frobenius_ratios_gaussian = zeros(length(mu_values), length(n_values));
frobenius_ratios_rademacher = zeros(length(mu_values), length(n_values));
thresholds_gaussian = zeros(length(mu_values), length(n_values));
thresholds_rademacher = zeros(length(mu_values), length(n_values));


for idx = 1:length(n_values)
    n = n_values(idx);
    
    %% Gaussian Matrix A
    A_gaussian = randn(n, p);
    columnNorms = vecnorm(A_gaussian, 2, 1).^2 / n;
    nu_gaussian = compute_rho(A_gaussian);
    
    %% Rademacher Matrix A
    A_rademacher = 2 * randi([0, 1], n, p) - 1;
    nu_rademacher = compute_rho(A_rademacher);
    
    
    for mu_idx = 1:length(mu_values)
        mu = mu_values(mu_idx);
        
        % Compute thresholds for each mu
        thresholds_gaussian(mu_idx, idx) = nu_gaussian / (1 + nu_gaussian);
        thresholds_rademacher(mu_idx, idx) = nu_rademacher / (1 + nu_rademacher);
             
        %% CVX Optimization for Gaussian
        W_e_gaussian=compute_W_e(A_gaussian);
        W_o_gaussian=compute_W_o(A_gaussian);
        
        % Frobenius norm ratio for Gaussian
        frobenius_ratios_gaussian(mu_idx, idx) = norm(W_e_gaussian - W_o_gaussian, 'fro') / norm(W_o_gaussian, 'fro');

        %% CVX Optimization for Rademacher
        W_e_rademacher=compute_W_e(A_rademacher);
        W_o_gaussian=compute_W_o(A_rademacher);
         
        % Frobenius norm ratio for Rademacher
        frobenius_ratios_rademacher(mu_idx, idx) = norm(W_e_rademacher - W_o_rademacher, 'fro') / norm(W_o_rademacher, 'fro');
    end
    
    % Display separate tables for each n and matrix type
    disp(['Results for Gaussian A, n = ', num2str(n)]);
    T_gaussian = table(mu_values', thresholds_gaussian(:, idx), frobenius_ratios_gaussian(:, idx), 'VariableNames', {'mu', 'nu_L_nu', 'Frobenius_Relative'});
    disp(T_gaussian);
    
    disp(['Results for Rademacher A, n = ', num2str(n)]);
    T_rademacher = table(mu_values', thresholds_rademacher(:, idx), frobenius_ratios_rademacher(:, idx), 'VariableNames', {'mu', 'nu_L_nu', 'Frobenius_Relative'});
    disp(T_rademacher);

end
