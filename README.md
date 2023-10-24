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


