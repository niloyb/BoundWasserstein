# BoundWasserstein

These scripts reproduce the results of the article **[Bounding Wasserstein distance with couplings](https://arxiv.org/abs/2112.03152)**  by Niloy Biswas and Lester Mackey. 

```
@article{biswas2021bounding,
  title={Bounding Wasserstein distance with couplings},
  author={Niloy Biswas and Lester Mackey},
  journal={arXiv preprint arXiv:2112.03152},
  year={2021}
}
```

All commands below should be run from this repository base directory.
- `estimators.R` contains functions for computing Wasserstein upper and lower bound estimates
    - `wp_ub_estimate(coupled_chain_sampler, no_chains, p=1, metric, parallel=TRUE)` returns coupled upper bound (CUB) estimates for the p-Wasserstein distance
    - `w2l2_lb_gelbrich(chain1,chain2)` returns the Gelbrich lower bound for 2-Wasserstein distance with the L2 metric
    - `W2LBmarginals(chain1,chain2)` returns the marginal-based lower bound for 2-Wasserstein distance with the L2 metric
    - `W2L2_UBLB_estimates(chain_sampler, no_chains, parallel=TRUE, lb='max_gelbrich_marginals')` returns 2-Wasserstein upper and lower bound estimates
- `kernels.R` contains functions for running Markov chains
    - `ula_proposal_mean(x, grad_log_pdf, sigma_mh)` returns the mean of the unadjusted Langevin algorithm (ULA) proposal given the current state x
    - `kernel_ula(current_state, grad_log_pdf,sigma_mh)` returns the next state in a ULA chain
    - `ula(init_state, grad_log_pdf, sigma_mh, iterations)` returns an ULA chain of length iterations
    - `kernel_mala(current_state, log_pdf, grad_log_pdf, sigma_mh)` returns the next state in a Metropolis-adjusted Langevin algorithm (MALA) chain
    - `mala(init_state, log_pdf, grad_log_pdf, sigma_mh, iterations)` returns a MALA chain of length iterations
    - `reflection_couping_normal_Idcov(mu1, mu2, norm_crn=NULL, unif_crn=NULL)` returns a reflection coupling of direct samplers for N(mu1, Id) and N(mu2, Id)
    - `reflection_couping_normal(mu1, mu2, cov_mat1_cholesky, cov_mat2_cholesky, cov_mat1_cholesky_inv, cov_mat2_cholesky_inv, norm_crn=NULL, unif_crn=NULL)` returns a reflection coupling of direct samplers for N(mu1, Sigma1) and N(mu2, Sigma2)
    - `coupled_kernel_ula(current_pi, current_q, pi_grad_log_pdf, q_grad_log_pdf, sigma_mh, reflect_threshold=0, norm_crn=NULL, unif_crn=NULL)` returns a coupling of ULA kernels targeting pi and q
    - `coupled_ula(init_pi, init_q, pi_grad_log_pdf, q_grad_log_pdf, sigma_mh, iterations, reflect_threshold=0)` returns a coupled ULA chain of length iterations
    - `coupled_kernel_mala(current_pi, current_q, pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf, sigma_mh, reflect_threshold=0, pi_correction=TRUE, q_correction=TRUE, norm_crn=NULL, unif_crn1=NULL, unif_crn2=NULL)` returns a coupling of MALA kernels targeting pi and q
    -  `coupled_mala(init_pi, init_q, pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf, sigma_mh, iterations, reflect_threshold=0, pi_correction=TRUE, q_correction=TRUE)` returns a coupled MALA chain of length iterations
    -  `leapfrog(x, v, gradlogtarget, stepsize, nsteps)` returns a leapfrog integrator update of position x and velocity v
    -  `hmc_kernel(current_state, current_pdf, logtarget, gradlogtarget, stepsize, nsteps)` returns the next state in a Hamiltonian Monte Carlo (HMC) chain
    -  `hmc(init_state, init_pdf, logtarget, gradlogtarget, stepsize, nsteps, chain_length)` returns an HMC chain of length chain_length
    -  `coupled_hmc_kernel(current_pi, current_q, current_pdf_pi, current_pdf_q, pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf, pi_stepsize, q_stepsize, pi_nsteps, q_nsteps, reflect_threshold, pi_correction, q_correction)` returns a coupling of HMC kernels targeting pi and q
    -  `coupled_hmc(init_pi, init_q, init_pdf_pi, init_pdf_q, 
           pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf,
           pi_stepsize, q_stepsize, pi_nsteps, q_nsteps,
           chain_length, reflect_threshold=0,
           pi_correction=TRUE, q_correction=TRUE)` returns a coupled MALA chain of length chain_length
-	`stylized_example_mvn` contains code for Figure 1.
    - To reproduce Figure 1 run `source('stylized_example_mvn/mvn_plots.R')`
-	`stylized_example_ula_mala` contains code for Figure 2.
    - To reproduce Figure 2 run `source('stylized_example_ula_mala/ula_mala_plots.R')`
-	`implementation_plots` contains code for Figure 3.
    - To reproduce Figures 3a and 3b run `source('implementation_plots/bimodal_plots.R')`
    - To reproduce Figure 3c run `source('implementation_plots/crn_reflection_plot.R')`
-	`huggins_etal_comparison` contains code for Figure 4 (left).
    - To reproduce Figure 4 (left) run `source('huggins_etal_comparison/viabel_simulations.R')`
-	`dobson_etal_comparison` contains code for Figure 4 (right).
    - To reproduce Figure 4 (right) run `source('dobson_etal_comparison/dobson_simulations.R')`
-	`bayesian_logistic_regression` contains code for Figure 5.
    - To reproduce Figure 5
        - Run `source('bayesian_logistic_regression/logreg_simulations.R')` with `data <- 'pima'` on [this line](https://github.com/niloyb/BoundWasserstein/blob/5801979f2022334ed4f6a4f9c4dd1a923bc93ce6/bayesian_logistic_regression/logreg_simulations.R#L41) to generate the Pima bounds dataframe `bayesian_logistic_regression/logreg_bounds_df_pima.RData`
        - Run `source('bayesian_logistic_regression/logreg_simulations.R')` with `data <- 'ds1'` on [this line](https://github.com/niloyb/BoundWasserstein/blob/5801979f2022334ed4f6a4f9c4dd1a923bc93ce6/bayesian_logistic_regression/logreg_simulations.R#L41) to generate the DS1 bounds dataframe `bayesian_logistic_regression/logreg_bounds_df_ds1.RData`
        - Run `source('bayesian_logistic_regression/logreg_plots.R')` to generate the figure plots
-	`half_t` contains code for Figure 6.
    - To reproduce Figure 6
        - Run `source('half_t/half_t_simulations.R')` to generate the bounds dataframes `half_t/half_t_bounds_df_riboflavin.RData` and `half_t/half_t_bounds_df_synthetic.RData`
        - Run `source('half_t/half_t_plots.R')` to generate the figure plots
-	`skinny_gibbs` contains code for Figure 7.
    - To reproduce Figure 7
        - Run `source('skinny_gibbs/skinny_gibbs_simulations.R')` to generate the bounds dataframe `skinny_gibbs/skinny_gibbs_df.RData`
        - Run `source('skinny_gibbs/skinny_gibbs_plots.R')` to generate the figure plot
-	`sinkhorn_comparison` contains code for supplementary Figure S1.
    - To reproduce supplementary Figure S1 run `source('sinkhorn_comparison/sinkhorn_simulations.R')`


