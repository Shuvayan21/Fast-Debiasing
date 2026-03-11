function [sens,spec]= W_e_sens_spec(y,A,sigma,beta,W_e)
    TG=zeros(p,1);
    tp=0;
    tn=0;
    fp=0;
    fn=0;
    infec=(abs(beta)>0);
    lambda = 1;
    cvx_begin quiet
                variable x_l(p)
                minimise (0.5*pow_pos(norm(y-A*x_l),2)+lambda*norm(x_l))
    cvx_end
    %W_e=compute_W_e(A);
    Sigma_beta_W=sigma^2/n*(W_e'*W_e);
    %taking out the delta part of the debiased estimate of x_d
    beta_l=x_l;
    beta_d_W=beta_l+1/(n)*W_e'*(y-A*beta_l);
    for i=1:1:p
        TG(i)=((beta_d_W(i)))/sqrt(Sigma_beta_W(i,i));
    end 
    J=(abs(TG)>2.33); %1% level of significance
    %sensitivity and specificity for the hypothesis test
    for i=1:1:n
        if(J(i)==1 && infec(i)==1)
            tp=tp+1;
        elseif (J(i)==1 && infec(i)==0)
            fp=fp+1;        
        elseif (J(i)==0 && infec(i)==0) 
            tn=tn+1;
        else
            fn=fn+1;
        end
    end
    sens=tp/(tp+fn);
    spec=tn/(tn+fp);
end