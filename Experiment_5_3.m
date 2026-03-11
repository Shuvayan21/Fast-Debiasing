rng(1)
n_values = 200:50:500; % n = 40, 50, 60, 70, 80, 90, 100
p = 500; % Number of columns
iterations = 1; % Set iterations to 1
 mu_values = linspace(0.2, 0.6, 41); % Various values of mu
 s=10;
 f_sig=0.05;
S = randperm(p, s); %support of the non-zero elements of \beta
    % Choosing non-zero elements from uniform distribution
    beta(S(1:floor(0.4*s))) =500+500*rand(floor(0.4*s),1);
    beta(S(floor(0.4*s)+1:s)) =500+500*rand(s-floor(0.4*s),1);
    %creating the standard deviation sigma
 frobenius_ratios_gaussian = zeros(length(mu_values), length(n_values));   
 sensitivity_W_o = zeros(length(mu_values), length(n_values)); 
 specificity_W_o = zeros(length(mu_values), length(n_values)); 
 sensitivity_W_e = zeros(length(mu_values), length(n_values)); 
 specificity_W_e = zeros(length(mu_values), length(n_values)); 
 elapsed_time_W_o = zeros(length(mu_values), length(n_values)); 
 elapsed_time_W_e = zeros(length(mu_values), length(n_values));
for idx = 1:length(n_values)
    n = n_values(idx);
    
    %% Gaussian Matrix A
    A = randn(n, p);
    rho = compute_rho(A);
    sigma=mean(abs(A*beta))*f_sig;
    eta = random("normal",0,sigma,[n 1]);
    y = A*beta+eta;
    
    for mu_idx = 1:length(mu_values)
        mu = mu_values(mu_idx);
        
             
        %% CVX Optimization for Gaussian
        tic
        W_e_gaussian=compute_W_e(A);
        elapsed_time_W_e(mu_idx,idx) = toc;
        tic
        W_o_gaussian=compute_W_o(A);
        elapsed_time_W_o(mu_idx,idx) = toc;
        % Frobenius norm ratio for Gaussian
        frobenius_ratios_gaussian(mu_idx, idx) = norm(W_e_gaussian - W_o_gaussian, 'fro') / norm(W_o_gaussian, 'fro');

        [sensitivity_W_e(mu_idx, idx),specificity_W_e(mu_idx, idx)] = W_e_sens_spec(y,A,sigma,beta,W_e_gaussian);
        [sensitivity_W_o(mu_idx, idx),specificity_W_o(mu_idx, idx)] = W_e_sens_spec(y,A,sigma,beta,W_o_gaussian);

    end
    

end
% Display separate tables for each n and mu
    disp(['Results for varying mu, n = ', num2str(n_values(5))]);
    T_mu = table(mu_values', sensitivity_W_e(:, 5), specificity_W_e(:, 5), sensitivity_W_o(:, 5), specificity_W_o(:, 5), elapsed_time_W_e(:,5), elapsed_time_W_o(:,5), frobenius_ratios_gaussian(:, 5), 'VariableNames', {'mu', 'sens_W_e', 'spec_W_e', 'sens_W_o', 'spec_W_o', 'time- W_e', 'time- W_o', 'Frobenius_Relative'});
    disp(T_mu);
%% Code for varying n  
frobenius_ratios_n = zeros( length(n_values),1);   
 sensitivity_W_o_n = zeros( length(n_values),1); 
 specificity_W_o_n = zeros( length(n_values),1); 
 sensitivity_W_e_n = zeros( length(n_values),1); 
 specificity_W_e_n = zeros( length(n_values),1); 
 elapsed_time_W_o_n = zeros( length(n_values),1); 
 elapsed_time_W_e_n = zeros( length(n_values),1); 
mu=zeros(n_values,1);
for idx = 1:length(n_values)
    n = n_values(idx);
    
    %% Gaussian Matrix A
    A = randn(n, p);
    rho = compute_rho(A);
    sigma=mean(abs(A*beta))*f_sig;
    eta = random("normal",0,sigma,[n 1]);
    y = A*beta+eta;
    
        mu(idx) = rho/(1+rho);
        
             
        %% CVX Optimization for Gaussian
        tic
        W_e_gaussian=compute_W_e(A);
        elapsed_time_W_e_n(idx) = toc;
        tic
        W_o_gaussian=compute_W_o(A);
        elapsed_time_W_o_n(idx) = toc;
        % Frobenius norm ratio for Gaussian
        frobenius_ratios_n(idx) = norm(W_e_gaussian - W_o_gaussian, 'fro') / norm(W_o_gaussian, 'fro');

        [sensitivity_W_e_n(idx),specificity_W_e_n(idx)] = W_e_sens_spec(y,A,sigma,beta,W_e_gaussian);
        [sensitivity_W_o_n(idx),specificity_W_o_n(idx)] = W_e_sens_spec(y,A,sigma,beta,W_o_gaussian);


    

end



    disp('Results for varying n ');
    T_n = table(n_values', sensitivity_W_e_n, specificity_W_e_n, sensitivity_W_o_n, specificity_W_o_n, elapsed_time_W_e_n, elapsed_time_W_o_n, frobenius_ratios__n, 'VariableNames', {'mu', 'sens_W_e', 'spec_W_e', 'sens_W_o', 'spec_W_o', 'time- W_e', 'time- W_o', 'Frobenius_Relative'});
    disp(T_n);