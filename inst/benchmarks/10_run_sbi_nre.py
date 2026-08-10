#!/usr/bin/env python
"""Train Python sbi's NRE on the shared simulations and sample each observation set.

Usage: python 10_run_sbi_nre.py --classifier resnet --n_samples 5000
Requires: pip install sbi (>= 0.22). sbi depends on nflows, which does not
build against setuptools >= 66 -- see README for the pin.

Both implementations get the same simulations, the same classifier family and
the same MCMC family (a slice sampler under build_posterior()'s defaults), so
the only thing that differs is the code under test. sbi's defaults are left
alone on purpose -- the point of the comparison is what a user gets out of the
box, not what either package can be tuned into.
"""
import argparse
import os

import numpy as np
import torch
from sbi.inference import NRE


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classifier", default="resnet",
                     choices=["resnet", "mlp", "linear"])
    ap.add_argument("--n_samples", type=int, default=5000)
    ap.add_argument("--prior_var", type=float, default=0.1)
    ap.add_argument("--noise_var", type=float, default=0.1)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    data_dir = os.path.join("data", "gaussian_linear_iid_nre")
    meta = [int(v) for v in open(os.path.join(data_dir, "meta.txt")).read().split()]
    dim, grid_n, n_obs = meta[0], meta[-1], meta[1:-1]

    theta = torch.tensor(np.loadtxt(os.path.join(data_dir, "theta.csv"),
                                    delimiter=","), dtype=torch.float32)
    x = torch.tensor(np.loadtxt(os.path.join(data_dir, "x.csv"),
                                delimiter=","), dtype=torch.float32)

    # task_gaussian_linear(): theta ~ N(0, prior_var I), x | theta ~ N(theta, noise_var I).
    prior = torch.distributions.MultivariateNormal(
        torch.zeros(dim), args.prior_var * torch.eye(dim)
    )

    inference = NRE(prior=prior, classifier=args.classifier)
    inference.append_simulations(theta, x).train()
    posterior = inference.build_posterior()

    out_dir = os.path.join("results", "gaussian_linear_iid_nre")
    os.makedirs(out_dir, exist_ok=True)
    for k in n_obs:
        x_obs = np.atleast_2d(np.loadtxt(
            os.path.join(data_dir, f"x_obs_n{k}.csv"), delimiter=","))
        samples = posterior.sample(
            (args.n_samples,), x=torch.tensor(x_obs, dtype=torch.float32)
        ).numpy()
        path = os.path.join(out_dir, f"sbi_{args.classifier}_n{k}.csv")
        np.savetxt(path, samples, delimiter=",")
        print("wrote", path)

    # Secondary metric (#146): the learned ratio itself, on the grid
    # 09_generate_data_nre.R wrote. potential_fn(theta) for an NRE-based
    # posterior is prior.log_prob(theta) + sum(log_ratio(theta, x_obs)), so
    # subtracting the (known, closed-form) prior term recovers the ratio up
    # to the same additive constant every NRE fit carries. This depends on
    # sbi's potential_fn/set_x API, which has moved between releases, so it
    # is wrapped defensively: a failure here should not cost the posterior
    # samples the C2ST comparison actually needs.
    grid_path = os.path.join(data_dir, "grid_theta.csv")
    x_obs1_path = os.path.join(data_dir, "x_obs_n1.csv")
    if grid_n > 0 and os.path.exists(grid_path) and os.path.exists(x_obs1_path):
        try:
            grid = torch.tensor(np.loadtxt(grid_path, delimiter=","),
                                dtype=torch.float32)
            if grid.ndim == 1:
                grid = grid.unsqueeze(1)
            x_obs1 = np.atleast_2d(np.loadtxt(x_obs1_path, delimiter=","))
            potential_fn = posterior.potential_fn
            potential_fn.set_x(torch.tensor(x_obs1, dtype=torch.float32))
            log_post = potential_fn(grid).detach()
            log_prior = prior.log_prob(grid)
            grid_ratio = (log_post - log_prior).numpy()
            path = os.path.join(out_dir, f"grid_log_ratio_sbi_{args.classifier}.csv")
            np.savetxt(path, grid_ratio, delimiter=",")
            print("wrote", path)
        except Exception as exc:  # depends on sbi's internal API; see comment above
            print(f"skipped grid log-ratio extraction: {exc}")


if __name__ == "__main__":
    main()
