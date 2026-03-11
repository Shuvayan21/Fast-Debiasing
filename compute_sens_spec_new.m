function [Sensitivity, Specificity] = compute_sens_spec_new(theta_true, D)
    % theta_true : vector of true coefficients (can be numeric)
    % D          : binary detection vector (1 = detected nonzero, 0 = detected zero)
    
    % Convert true coefficients to binary (1 = nonzero, 0 = zero)
    true_binary = theta_true ~= 0;
    
    % True positives, false negatives
    TP = sum((true_binary == 1) & (D == 1));
    FN = sum((true_binary == 1) & (D == 0));
    
    % True negatives, false positives
    TN = sum((true_binary == 0) & (D == 0));
    FP = sum((true_binary == 0) & (D == 1));
    
    % Compute Sensitivity and Specificity
    Sensitivity = TP / (TP + FN);
    Specificity = TN / (TN + FP);
end
