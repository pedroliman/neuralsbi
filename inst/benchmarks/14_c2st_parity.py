#!/usr/bin/env python
"""Score the sample-set pairs 13_c2st_parity.R wrote with sbibm's own C2ST.

Usage: python 14_c2st_parity.py
Requires: pip install scikit-learn (sbibm itself is not needed -- the function
below is sbibm/metrics/c2st.py with the torch tensor calls written in numpy,
which is the same arithmetic on the same defaults and drops a heavy dependency
from a check that does not use a network).

Prints the two tables side by side. The two implementations initialize and
shuffle from different random number streams, so the columns agree to
Monte-Carlo noise rather than to the digit; a gap of more than about 0.02 is
worth looking into.
"""
import csv
import os

import numpy as np
from sklearn.model_selection import KFold, cross_val_score
from sklearn.neural_network import MLPClassifier

CASES = ["same_d2", "shift_d2", "scale_d2", "shift_d1", "scale_d5"]
RES_DIR = os.path.join("results", "c2st_parity")


def c2st(X, Y, seed=1, n_folds=5, scoring="accuracy", z_score=True):
    if z_score:
        X_mean = X.mean(axis=0)
        X_std = X.std(axis=0, ddof=1)  # torch.std is unbiased by default
        X = (X - X_mean) / X_std
        Y = (Y - X_mean) / X_std
    ndim = X.shape[1]
    clf = MLPClassifier(
        activation="relu",
        hidden_layer_sizes=(10 * ndim, 10 * ndim),
        max_iter=10000,
        solver="adam",
        random_state=seed,
    )
    data = np.concatenate((X, Y))
    target = np.concatenate((np.zeros(X.shape[0]), np.ones(Y.shape[0])))
    shuffle = KFold(n_splits=n_folds, shuffle=True, random_state=seed)
    return float(np.mean(cross_val_score(clf, data, target, cv=shuffle, scoring=scoring)))


def read(name, side):
    path = os.path.join(RES_DIR, f"{name}_{side}.csv")
    return np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)


def main():
    ours = {}
    ours_path = os.path.join(RES_DIR, "neuralsbi.csv")
    if os.path.exists(ours_path):
        with open(ours_path, newline="") as fh:
            for row in csv.DictReader(fh):
                ours[row["case"]] = (float(row["accuracy"]), float(row["auc"]))
    else:
        print(f"{ours_path} not found -- run 13_c2st_parity.R first")

    print(f"{'case':10s} {'sbibm acc':>10s} {'ours acc':>9s} "
          f"{'sbibm auc':>10s} {'ours auc':>9s}")
    for name in CASES:
        X, Y = read(name, "x"), read(name, "y")
        acc = c2st(X, Y, seed=1)
        auc = c2st(X, Y, seed=1, scoring="roc_auc")
        r_acc, r_auc = ours.get(name, (float("nan"), float("nan")))
        print(f"{name:10s} {acc:10.4f} {r_acc:9.4f} {auc:10.4f} {r_auc:9.4f}")


if __name__ == "__main__":
    main()
