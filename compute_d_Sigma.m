function W= compute_d_Sigma(A,Sig,beta)
    [n,p]=size(A);
    d=1/(1-nnz(beta)/p);
    W=d*inv(Sig)*A';
end