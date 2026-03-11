# Fast-Debiasing
Matlab functions for the implementation of Fast Debiasing and obtaining the simulation results of the paper "Fast Debiasing of the LASSO Estimator" accepted in TMLR.

Description of Codes:
compute_rho.m - Given a matrix A, computes the quatity rho given in Theorem 1.
compute_W_e.m - Computes the exact solution W_e provided in the Fast Debiasing work given a design matrix A.
compute_W_o.m - Computes the optimal solution W_o by solving the optimization problem given in Eqn.(6) of the paper given a design matrix A.
compute_d_Sigma - Computes the inverse M=d\Sigma^-1 for Sec 5.2 of the paper.
W_e_sens_spec.m - Given measurements y, design matrix A and signal \beta, performs debiasing using W_e and returns the sensitivity and specificity of the debiased lasso estimates.
W_o_sens_spec.m - Given measurements y, design matrix A and signal \beta, performs debiasing using W_o and returns the sensitivity and specificity of the debiased lasso estimates.
d_Sig_sens_spec.m - Given measurements y, design matrix A and signal \beta, performs debiasing using M=D\Sigma^-1 and returns the sensitivity and specificity of the debiased lasso estimates.
compute_sens_spec_new.m - Given the ground truth \beta and an estimated signal \hat{\beta}, returns the sensitivity and specificity by maximizing the Youden's index.
generate_crosstalk_matrix.m- Creates a Rademacher matrix with crosstalk/correlated rows as described in Sec~5.2.
makeSigma.m - Creates a banded matrix as described in Sec. 5.2.
plot_colored_cube_single_figure.m - Code to plot hyperspectral images.
Experiment_5_<a>.m - Generates results for Sec 5.<a> where a=1,2,3,4.
Experiment_6_1_uncorrelated.m - Generates results for Sec. 6.1 for uncorrelated entries in the design matrix A.
Experiment_6_1_correlated.m - Generates results for Sec. 6.1 for crosstalk in the design matrix A.
Experiment_6_2_extraction.m - Extracts 30 sets of 5 frames each from Waterfall.mp4 video.
Experiment_6_2_reconstruction.m - Performs multi-frame video reconstruction using LASSO, Debiased LASSO on the extracted frames of waterfall video.
Experiment_6_2_video_creation.m - Creates a comparison video based on the multi-frame video reconstruction.
Experiment_6_3.m - Generates results for Sec 6.3.

Experiment_6_2_data_and_results - Contains the original waterfall.mp4 and the comparison video for different reconstructions using LASSO and Debiased LASSO along with the extracted and reconstructed frames.
