#!/usr/bin/env python
"""Train Python sbi's NPE on the shared simulations and save posterior samples.

Usage: python 02_run_sbi_python.py --task gaussian_linear --estimator maf
Requires: pip install sbi pandas (sbi >= 0.22)
"""
import argparse
import os

import numpy as np
import torch
from sbi.inference import NPE
from sbi.utils import BoxUniform

PRIORS = {
    "gaussian_linear": lambda: torch.distributions.MultivariateNormal(
        torch.zeros(10), 0.1 * torch.eye(10)
    ),
    "two_moons": lambda: BoxUniform(-torch.ones(2), torch.ones(2)),
    "slcp": lambda: BoxUniform(-3 * torch.ones(5), 3 * torch.ones(5)),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", default="gaussian_linear")
    ap.add_argument("--estimator", default="maf", choices=["maf", "mdn", "nsf"])
    ap.add_argument("--n_samples", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    torch.manual_seed(args.seed)
    data_dir = os.path.join("data", args.task)
    theta = torch.tensor(np.loadtxt(os.path.join(data_dir, "theta.csv"),
                                    delimiter=","), dtype=torch.float32)
    x = torch.tensor(np.loadtxt(os.path.join(data_dir, "x.csv"),
                                delimiter=","), dtype=torch.float32)
    x_obs = np.atleast_2d(np.loadtxt(os.path.join(data_dir, "x_obs.csv"),
                                     delimiter=","))

    inference = NPE(prior=PRIORS[args.task](), density_estimator=args.estimator)
    inference.append_simulations(theta, x).train()
    posterior = inference.build_posterior()

    out_dir = os.path.join("results", args.task)
    os.makedirs(out_dir, exist_ok=True)
    for i, xo in enumerate(x_obs, start=1):
        samples = posterior.sample(
            (args.n_samples,), x=torch.tensor(xo, dtype=torch.float32)
        ).numpy()
        path = os.path.join(out_dir, f"sbi_{args.estimator}_obs{i}.csv")
        np.savetxt(path, samples, delimiter=",")
        print("wrote", path)


if __name__ == "__main__":
    main()
