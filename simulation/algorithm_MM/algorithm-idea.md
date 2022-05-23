---
output:
  html_document: default
  pdf_document: default
---
## Introduction

Consider a scalar random variable $Y$ and a $D$-dimensional random vector $X$ such that the following linear quantile regression model holds
$$
Q(\tau | x) = x'\beta(\tau),
$$
for $\tau\in(0,1)$ and $x\in\mathrm{support}(X)$, where $\beta:(0,1)\to\mathbb{R}^D$ is a smooth function and where
$$
Q(\tau|x):=\inf\{y\in\mathbb{R} : F(y|x)\ge\tau\},
$$
with $F(\cdot|x)$ denoting the conditional cdf of $Y$ given $X=x$. It is well known that, given any vector $b\in\mathbb{R}^D$ and $\tau\in(0,1)$, one has
$$
\mathbb{E} \rho_\tau \big(Y - X'b\big) \ge \mathbb{E} \rho_\tau \big(Y - X'\beta(\tau)\big),
$$
where $2\rho_\tau(z) := (2\tau-1)z + |z|$ for $z\in\mathbb{R}$ and $\tau\in(0,1)$. Hence, given a grid
$$
0<\tau_1\le\tau_2\le\cdots\le\tau_M<1
$$
of quantile levels (possibly depending on $N\in \mathbb{N}$), and writing
$$
\boldsymbol{\beta} = \begin{pmatrix}
\beta_1(\tau_1) & \beta_1(\tau_2) & \cdots & \beta_1(\tau_M)\\
\beta_2(\tau_1) & \beta_2(\tau_2) & \cdots & \beta_2(\tau_M)\\
\vdots & \vdots & \ddots & \vdots\\
\beta_D(\tau_1) & \beta_D(\tau_2) & \cdots & \beta_D(\tau_M)\\
\end{pmatrix},
$$
it also holds that, for any $D\times M$ matrix $B$ with columns $b_m\in\mathbb{R}^D$ ($m=1,\dots,M$),
$$
R(B):= \sum_{m=1}^M\mathbb{E} \rho_{\tau_m}\big(Y - X'b_m\big) \ge \sum_{m=1}^M\mathbb{E} \rho_{\tau_m}\big(Y - X'\beta(\tau_m)\big) = R(\boldsymbol{\beta})
$$
Said another way, we have
$$
R(\boldsymbol{\beta}) = \arg\min_B R(B)
$$
where the minimum in $B$ runs through all $D\times M$ matrices $B$.

Now let $\varphi_1,\varphi_2,\dots,\varphi_L$ be linearly independent vectors in $\mathbb{R}^M$ (possibly depending on $N\in\mathbb{N}$), with $L\le M$. Also assume that
$$
\boldsymbol{\beta}\approx \boldsymbol{\alpha\varphi}\qquad\qquad\tag{proj}
$$
for some $D\times L$ matrix $\boldsymbol\alpha$,  where $\boldsymbol{\varphi}$ is the $L\times M$ matrix with _rows_ $\varphi_\ell$. This will be the case whenever the mapping $\tau\mapsto \beta(\tau)$ is "sufficiently" regular, in the sense that, for all $d$, the coefficient $\beta_d(\cdot)$ can be represented in terms of a sequence of basis functions $f_1, f_2, \dots$ as
$$
\beta_d(\cdot) = \sum_{\ell=1}^\infty \alpha_{d\ell}f_\ell(\cdot).
$$

In this setting, one could take $[\varphi_{\ell}]_m = f_\ell(\tau_m)$ for example. To be concrete, the vectors $\varphi_1,\dots,\varphi_L$ can be taken from application of a Gram-Schmidt algorithm to the set $\{\tilde{\varphi}_1,\dots,\tilde{\varphi}_L\}\subseteq\mathbb{R}^M$ given by
$$
\tilde{\varphi}_\ell' = \begin{pmatrix}
\tau_1^{\ell-1} & \cdots & \tau_M^{\ell-1}
\end{pmatrix},\qquad \qquad\ell=1,\dots,L.
$$

## Estimation

Given a random sample $(X_n,Y_n)_{n=1}^N$ from $(X,Y)$, the representation in equation $(\textrm{proj})$ suggests estimating $\boldsymbol{\beta}$ by finding a $D\times L$ matrix $\widehat{\boldsymbol\alpha}$ to solve the optimization problem
$$
\min_{\boldsymbol{a}} \sum_{m=1}^M\sum_{n=1}^N\rho_{\tau_m}(Y_n - X_n'\boldsymbol{a}\varphi_m)
$$
where $\varphi_m$ is the $m$th *column* of $\boldsymbol{\varphi}$. Letting $\boldsymbol{\eta}$ denote an $N\times M$ matrix with positive entries, we can use the [eta-trick](https://francisbach.com/the-η-trick-or-the-effectiveness-of-reweighted-least-squares/) to write the above minimization problem alternatively as
$$
\min_{\boldsymbol {a}} \min_{\boldsymbol{\eta}>0} \widehat{R}(\boldsymbol{a},\boldsymbol{\eta})
$$
where
$$
2\widehat{R}(\boldsymbol{a},\boldsymbol{\eta}) := \sum_{m=1}^M\sum_{n=1}^N (2\tau_m-1)(Y_n-X_n'\boldsymbol{a}\varphi_m) + \frac12\left[\frac{(Y_n - X_n'\boldsymbol{a}\varphi_m)^2}{\eta_{nm}} + \eta_{nm}\right]
\tag{1}
$$
and
$$
\eta_{nm} = |Y_n - X_n'\boldsymbol{a}\varphi_m|.\qquad\qquad\tag{2}
$$
This can be extended to include group lasso-type penalties, for example by letting $\boldsymbol{\zeta}$ denote a $D\times M$ matrix and writing
$$
\widehat{R}_\lambda(\boldsymbol{a},\boldsymbol{\eta},\boldsymbol{\zeta}) = \widehat{R}(\boldsymbol{a},\boldsymbol{\eta}) + \lambda\sum_{d=1}^D\frac12\left[\frac{\Vert \boldsymbol{a}_{d}\Vert ^2}{\zeta_{d}}+\zeta_{d}\right]
\tag{3}
$$
where $\lambda>0$ is a tuning parameter,
$$
\zeta_{d} = \Vert \boldsymbol{a}_{d} \Vert
$$
and
$$
\Vert \boldsymbol{a}_d \Vert =  \sqrt{\sum_{\ell = 1}^L{\boldsymbol{a}_{d\ell}^2}}
$$

Now, let 

$$
vec(\boldsymbol{a}) = [a_{11}, a_{21}, \dots, a_{d\ell}],
$$

$$
vec(X_n \varphi_m') = \begin{bmatrix}
                        X_{n1} \varphi_{1m} \\ 
                        X_{n2} \varphi_{1m} \\
                        \vdots \\
                        X_{nd} \varphi_{1m} \\
                        X_{n1} \varphi_{2m} \\
                        X_{n1} \varphi_{2m} \\
                        \vdots \\
                        X_{n1} \varphi_{\ell m} \\
                        X_{n2} \varphi_{\ell m} \\
                        \vdots \\
                        X_{nd} \varphi_{\ell m}
                        \end{bmatrix}
$$

With this, we re-write $X_n'\boldsymbol{a}\varphi_m$ in equation $(1)$ by $vec(X_n \varphi_m')'vec(\boldsymbol{a}) = vec(\boldsymbol{a})'vec(X_n \varphi_m')$.

In this notation, we differentiate $(3)$ and find the estimator for $\boldsymbol{a}_{vec}$ as follows: 

$$
\widehat{vec(\boldsymbol{a})} = \left[\sum_{n=1}^{N} \sum_{m=1}^{M} \left[- \frac{1}{\eta_{nm}}vec(X_n \varphi_m')vec(X_n \varphi_m')' + \lambda diag(vec(\frac{1}{\zeta})) \right]\right]^{-1} + \\
\left[\sum_{n=1}^{N} \sum_{m=1}^{M} \left[(1-\tau_m) vec(X_n \varphi_m')- \frac{Y_n}{\eta_{nm}}vec(X_n \varphi_m')\right]\right]
\tag{4}
$$

## Algorithm

0. Initialize $vec(\boldsymbol{a})^0$ arbitrarily.
1. Given the value $vec(\boldsymbol{a})^r$ of the $r$-th iteration, find a solution $vec(\boldsymbol{a})^{r+1}$ for equation $(4)$ with $\eta_{nm}^r$ in place of $\eta_{nm}$ and $\zeta^r$ in place of $\zeta$. 
2. Use a precision parameter $\pi$ to set $\eta_{nm}^{r} = \max\{\pi,\eta_{nm}^{r}\}$ and $\zeta^{r} = \max\{\pi,\zeta^{r}\}$.
3. Iterate steps 1-2 until the convergence criterion has been met: $\frac{abs(R_{\lambda}(vec(\boldsymbol{a})^{r},\eta^{r},\zeta^{r}) - R_{\lambda}(vec(\boldsymbol{a})^{r+1},\eta^{r+1},\zeta^{r+1}))}{R_{\lambda}(vec(\boldsymbol{a})^{r},\eta^{r},\zeta^{r})} < tol$.

