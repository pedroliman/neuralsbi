#!/usr/bin/env python
"""Train Python sbi's NLE on the shared simulations and sample each observation set.

Usage: python 06_run_sbi_nle.py --estimator maf --n_samples 5000
Requires: pip install sbi (>= 0.23)

Both implementations get the same simulations, the same estimator family and
the same MCMC family, so the only thing that differs is the code under test.
sbi's defaults are left alone on purpose -- the point of the comparison is what
a user gets out of the box, not what either package can be tuned into.
"""
import argparse
import os

import numpy as np
import torch
from sbi.inference import NLE


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimator", default="maf", choices=["maf", "mdn", "nsf"])
    ap.add_argument("--n_samples", type=int, default=5000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    data_dir = os.path.join("data", "gaussian_linear_iid")
    meta = [int(v) for v in open(os.path.join(data_dir, "meta.txt")).read().split()]
    dim, n_obs = meta[0], meta[1:]

    theta = torch.tensor(np.loadtxt(os.path.join(data_dir, "theta.csv"),
                                    delimiter=","), dtype=torch.float32)
    x = torch.tensor(np.loadtxt(os.path.join(data_dir, "x.csv"),
                                delimiter=","), dtype=torch.float32)

    # task_gaussian_linear(): theta ~ N(0, 0.1 I), x | theta ~ N(theta, 0.1 I).
    prior = torch.distributions.MultivariateNormal(
        torch.zeros(dim), 0.1 * torch.eye(dim)
    )

    inference = NLE(prior=prior, density_estimator=args.estimator)
    inference.append_simulations(theta, x).train()
    posterior = inference.build_posterior()

    out_dir = os.path.join("results", "gaussian_linear_iid")
    os.makedirs(out_dir, exist_ok=True)
    for k in n_obs:
        x_obs = np.atleast_2d(np.loadtxt(
            os.path.join(data_dir, f"x_obs_n{k}.csv"), delimiter=","))
        samples = posterior.sample(
            (args.n_samples,), x=torch.tensor(x_obs, dtype=torch.float32)
        ).numpy()
        path = os.path.join(out_dir, f"sbi_{args.estimator}_n{k}.csv")
        np.savetxt(path, samples, delimiter=",")
        print("wrote", path)


if __name__ == "__main__":
    main()
