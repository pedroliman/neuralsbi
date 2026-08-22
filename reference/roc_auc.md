# Area under the ROC curve

The rank form of the Mann-Whitney statistic, which needs no threshold
sweep and handles ties the way `scikit-learn`'s `roc_auc_score` does.

## Usage

``` r
roc_auc(prob, label)
```

## Arguments

- prob:

  Predicted probability of class 1.

- label:

  0/1 labels.
