# Polynomial Interpolation Strategies

## Overview

Implements and evaluates three polynomial interpolation forms — Barycentric-1,
Barycentric-2, and Newton's Form — across three mesh distributions (uniform,
Chebyshev of the first kind, Chebyshev of the second kind) and three node
orderings (increasing, decreasing, Leja). All stability and accuracy experiments
are run in single precision, with double precision serving as the reference value
for error analysis. Four test functions expose distinct aspects of conditioning,
numerical stability, and convergence behavior.

## Methodology

Each interpolating form represents the same polynomial $p_n(x)$ through a
distinct mathematical formulation. The forms are implemented as paired routines:
a coefficient generation step and a polynomial evaluation step.

**Barycentric-1** computes weights $\gamma_i = 1/\omega'_{n+1}(x_i)$ via a
nested product loop — $O(n^2)$ time, $O(n)$ space — then evaluates:

$$
p_n(x) = \omega_{n+1}(x) \sum_{i=0}^{n} \frac{y_i \gamma_i}{x - x_i}
$$

**Barycentric-2** exploits closed-form weight formulas for each mesh type,
reducing coefficient generation from $O(n^2)$ to $O(n)$. Evaluation is $O(n)$
per point via:

$$
p_n(x) = \frac{\displaystyle\sum_{i=0}^{n} \frac{\beta_i y_i}{x - x_i}}
               {\displaystyle\sum_{i=0}^{n} \frac{\beta_i}{x - x_i}}
$$

Weights for each mesh: uniform $\beta_i = (-1)^i \binom{n}{i}$, Chebyshev-1
$\beta_i = (-1)^i \sin\!\left(\frac{(2i+1)\pi}{2n+2}\right)$, Chebyshev-2
$\beta_i = (-1)^i \delta_i$ where $\delta_i = 1/2$ at endpoints and $1$
elsewhere.

**Newton's Form** generates divided difference coefficients via an in-place
backward update — $O(n^2)$ time, $O(n)$ space — then evaluates using
adapted Horner's rule in $O(n)$ operations:

$$
p_n(x) = \sum_{i=0}^{n} y[x_0,\ldots,x_i]\,\omega_i(x),
\qquad s \leftarrow s(x - x_i) + \alpha_i
$$

Node ordering for Newton's Form is controlled by a separate routine supporting
increasing, decreasing, and Leja ordering. Leja ordering selects nodes
sequentially to maximize $\prod_{j<k}|x - x_j|$, preventing basis function
overflow in single precision for high-degree or clustered meshes.

Conditioning is measured by two quantities computed via a dedicated routine
forced into double precision:

$$
\Lambda_n = \max_{x \in [a,b]} \kappa(x,n,1) = \max_x \sum_i |l_i(x)|,
\qquad
H_n = \max_{x \in [a,b]} \kappa(x,n,y) = \max_x \frac{|\omega_{n+1}(x)
\sum \gamma_i f_i / (x-x_i)|}{|p_n(x)|}
$$

Error in single precision is bounded by:

$$
\|f(x) - p_n(x)\|_\infty \leq H_n \mu + \|f - p_n\|_\infty, \qquad
\mu = \text{unit round-off}
$$

## Test Functions

**$f_1(x) = (x-2)^9$, $n=9$, $[-1,1]$:** Chebyshev meshes produce
significantly lower $\Lambda_n$ and $H_n$ than uniform. Barycentric-2 is
consistently stable across all meshes. Leja ordering reduces Newton's Form
error on uniform mesh relative to increasing/decreasing order.

**$f_2(x) = \prod_{i=1}^d (x-i)$, $[0, d+1]$, $d = 5, 20, 29$:** The
problem is intrinsically ill-conditioned with $H_n \sim 10^{16}$ regardless of
mesh. Increasing and decreasing Newton orderings fail with catastrophic overflow
at $d = 20$ and $d = 29$. Leja ordering and Barycentric-2 reduce relative error
by several orders of magnitude.

**$f_3(x) = \ell_n(x)$ (Lagrange basis), $n = 20, 29$, $[-1,1]$:**
$H_n = 1$ throughout since the function only takes values 0 and 1. $\Lambda_n$
grows rapidly on a uniform mesh due to the Runge phenomenon while remaining
near-constant on Chebyshev meshes. On a uniform mesh, increasing order is the
best Newton performer; Leja ordering deteriorates badly at $n = 29$.

**$f_4(x) = 1/(1+25x^2)$ (Runge function), $n = 5, \ldots, 100$, $[-1,1]$:**
The uniform mesh diverges — error reaches $\infty$ at $n = 50$ and produces
NaN for $n \geq 60$. Both Chebyshev meshes converge monotonically. The
precision threshold is reached near $n = 100$ at $\|r\|_\infty \approx
2.23 \times 10^{-7}$, approaching single precision unit round-off. Leja-ordered
Newton on Chebyshev nodes matches Barycentric-2 accuracy; increasing/decreasing
Newton on Chebyshev nodes fails at $n = 30$.

## Language

MATLAB

## How to Run

1. Only need `Polynomial_Driver.m` file
2. Run the main test driver — it executes all four test functions across all
   method, mesh, and ordering combinations
3. Each run outputs: infinity norm error, relative infinity norm, Lebesgue
   constant $\Lambda_n$, conditioning number $H_n$, RMS error, and node error
4. Interpolation plots and residual plots are generated for each configuration
5. The convergence plot for $f_4$ is produced separately showing
   $\|f - p_n\|_\infty$ vs. $n$ for all three mesh types alongside the
   single precision floor
