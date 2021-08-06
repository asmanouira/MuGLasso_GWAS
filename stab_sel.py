import numpy as np


def stability_selection(MugLasso, alphas, n_bootstrap_iterations,
                        X, y, seed):
  n_samples, n_variables = X.shape
  n_alphas = alphas.shape[0]

  rnd = np.random.RandomState(seed)
  selected_variables = np.zeros((n_variables,
                                 n_bootstrap_iterations))
  stability_scores = np.zeros((n_variables, n_alphas))

  for idx, alpha, in enumerate(alphas):
    # This is the sampling step, where bootstrap samples are generated
    
    for iteration in range(n_bootstrap_iterations):
      bootstrap1 = rnd.choice(np.arange(n_samples_pop1),
                             size=n_samples_pop1 // 2,
                             replace=False)
      bootstrap2 = rnd.choice(np.arange(n_samples_pop2),
                             size=n_samples_pop2 // 2,
                             replace=False)
      bootstrap = bootstrap1 + bootstrap2

      X_train = X[bootstrap, :]
      y_train = y[bootstrap]

      # Assume scikit-learn implementation
      MugLasso.set_params({'C': alpha}).fit(X_train, y_train)
      selected_variables[:, iteration] = (np.abs(MugLasso.coef_) > 1e-4)


    stability_scores[:, idx] = selected_variables.mean(axis=1)

  return stability_scores