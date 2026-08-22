# Cross-validated logistic regression for [`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md)

Cross-validated logistic regression for
[`c2st()`](https://neuralsbi.pedrodelima.com/reference/c2st.md)

## Usage

``` r
c2st_logistic_prob(x_train, y_train, x_test)
```

## Arguments

- x_train, y_train:

  Training draws and their 0/1 labels.

- x_test:

  Draws to score.

## Value

Predicted probability of class 1 for each row of `x_test`.
