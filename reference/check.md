# Argument validation shared across the package

One place for the checks that public entry points repeat. Every message
names the argument it is about, says what was actually wrong, and where
there is a relevant help topic it points at it. That is the voice of
[`as_sim_draw()`](https://neuralsbi.pedrodelima.com/reference/as_sim_draw.md),
which is the model these follow: a user who mistypes an argument should
not have to guess which of `theta` and `x` the complaint is about, or
read a coercion error raised three frames down.

## Details

These are internal helpers, not exported. Call them at public
boundaries; internal code that already knows its shapes keeps using
[`as_theta_matrix()`](https://neuralsbi.pedrodelima.com/reference/as_theta_matrix.md).
