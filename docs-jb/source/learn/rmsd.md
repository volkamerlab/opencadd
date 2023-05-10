# RMSD

Root-mean-square deviation (RMSD), also called root-mean-square error (RMSE),
is a statistical measure of difference between two instances of a set of datapoints.
For example, one set may be the predicted values of a quantity (e.g. obtained by a regression model), whereas the
other set represents the true observed values. In this case, the RMSD acts as a measure of accuracy 
for the prediction. Similarly, the dataset may be the positions of atoms in a molecule, where each
instance corresponds to one conformation. Here, RMSD can act as a measure of similarity between the
two conformations. 

In general, a dataset $\hat{\textbf{P}}$ can be defined as an $N\times K$ matrix, 
containing $N$ points $\hat{\textbf{p}}_i$ (rows) in a $K$-dimensional space:

```{math}
:label: P
\hat{\textbf{P}}_{N\times K} 
= \begin{bmatrix}
    \hat{\textbf{p}}_1 \\
    \vdots \\
    \hat{\textbf{p}}_{N}
    \end{bmatrix}
= \begin{bmatrix}
    \hat{p}_{11} & \dots & \hat{p}_{1K} \\
    \vdots & \ddots & \vdots \\
    \hat{p}_{N1} & \dots & \hat{p}_{NK}
    \end{bmatrix}
```

```{note}
For scalar values, $K$ is equal to 1, and $\hat{\textbf{P}}$ can be thought of as an $N\times 1$ matrix.
```

Given another instance of the same data as another $N\times K$ matrix $\textbf{P}$,
and assuming a one-to-one correlation between the points in $\hat{\textbf{P}}$ and $\textbf{P}$
(i.e. assuming the values in both datasets are in identical order),
the RMSD is a function that maps the two datasets to a single real value 
in the interval $\left[0,\infty\right)$:

```{math}
\textrm{RMSD}: \mathbb{R}^{N \times K} \times \mathbb{R}^{N \times K} \rightarrow \mathbb{R}_{\ge 0}
```

It is defined as:

```{math}
:label: RMSD
\textrm{RMSD}(\hat{\textbf{P}}, \textbf{P})
= \sqrt{\frac{1}{N}\sum_{n = 1}^{N}\|\hat{\textbf{p}}_n - \textbf{p}_n\|^2} \\
```

```{hint}
There is no differentiation between $\hat{\textbf{P}}$ and $\textbf{P}$, i.e.
$\textrm{RMSD}(\hat{\textbf{P}}, \textbf{P}) = \textrm{RMSD}(\textbf{P}, \hat{\textbf{P}})$
```

The term $\|\hat{\textbf{p}}_n - \textbf{p}_n\|$ in the above equation corresponds to the *D* in RMS**D**, 
i.e. the *deviation*, also known as *error* or *residual*. It is the euclidean distance (i.e. 2-norm) 
between two correlated datapoints:

```{math}
:label: 2-norm
\|\hat{\textbf{p}}_n - \textbf{p}_n\|
= \sqrt{\sum_{k = 1}^{K}\left(\hat{p}_{nk} - p_{nk}\right)^2} \\
```

Plugging eq. {eq}`2-norm` into {eq}`RMSD`, the RMSD function can also be rewritten as:

```{math}
:label: RMSD2
\textrm{RMSD}(\hat{\textbf{P}}, \textbf{P})
= \sqrt{\frac{1}{N}\sum_{n = 1}^{N}\sum_{k = 1}^{K}\left(\hat{p}_{nk} - p_{nk}\right)^2}
```

Put into words, the RMSD between two sets of one-to-one correlated ($K$-dimensional) points is
obtained by first calculating the pairwise *deviation*, i.e. the euclidean distance, between each
pair of points, squaring each distance, and calculating the arithmetic mean of the squared distances.
The RMSD value is then the square root of this arithmetic mean. In other words, it is the square root 
of the mean of squared deviations (distances). Therefore, when two datasets are completely identical,
all pairwise distances become zero, resulting in an RMSD value of zero. On the other hand, 
there is no upper limit for the RMSD value; as the distances between two datasets tends to infinity,
so does the RMSD.

```{important}
The main piece of information being calculated by the RMSD function is the pairwise euclidean distance of
correlated points in the two datasets; the rest is concerned with averaging all these pairwise distances 
to a single scalar metric representing the difference between the two sets. The euclidean distance is a 
sensible metric for representing difference, when what matters is the absolute positions of points in the datasets 
and not their relaive orientations to one another. This is the case, e.g. when the datasets represent the 
predicted and observed values of a function. On the other hand, consider the case where RMSD is being 
used to assess the structural difference between two objects, e.g. two conformations of the same molecule. 
Thus, each point in the dataset would correspond to the position of a specific point in the object,
for example, position of each atom in the molecule. In this case, simply measuring the euclidean distance of
pairs of correlated atoms will not necessarily provide a sensible measure for structural difference 
between the two conformations. The reason is that the similarity of two objects does not depend on their absolute
position in physical space, but on the relative orientation of characteristic points within those objects.
Imagine having two identical statues in front of you. No matter how you move each one of them in space,
they remain structurally identical, and your brain is capable to realize that as well. In a sense, the brain is
first superposing the two structures on top of each other as much as possible, and then assessing their
similarity based on how much they overlap. In contrast, the RMSD between two identical objects will increase
as they are simply moved away from each other in space, since we are measuring euclidean distances between the
points. It is only when all points of the two identical objects are perfectly superposed on each other, 
that the RMSD becomes zero, indicating correctly that the two objects are indeed structurally identical.
Therefore, in order to use RMSD to correctly measure the structural difference between two objects 
represented by a set of points, it is crucial to first align the two objects as good as possible. 
```


````{note}
RMSD is similar to other regression metrics, such as the mean absolute error (MAE), which is
defined as the arithmetic mean of the absolute values of pairwise distances:
```{math}
:label: MAE
\textrm{MAE}(\hat{\textbf{P}}, \textbf{P})
= \frac{1}{N}\sum_{n = 1}^{N}\lvert\hat{\textbf{p}}_n - \textbf{p}_n\rvert
```
However, whereas the MAE value is proportional to the distances (errors), RMSD is proportional 
to squared distances. Consequently, the RMSD value is more sensitive to outliers, 
as it quadratically penalizes the errors.  

Another metric is the mean square error (MSE), which lies in between RMSD and MAE in terms of similarity.
Instead of calculating the arithmetic mean of absolute errors (what MAE does), it calculates the
arithmetic mean of squared errors:
```{math}
:label: MSE
\textrm{MSE}(\hat{\textbf{P}}, \textbf{P})
= \frac{1}{N}\sum_{n = 1}^{N}\left(\hat{\textbf{p}}_n - \textbf{p}_n\right)^2
```
In case of MSE, the errors are also penalized quadratically as in RMSD. However, since the final 
square root is missing, the unit of MSE values are the squared units of the dataset values. 
This makes the interpretation of MSE values less straightforward, in constrast to the RMSD values, 
which are in the same unit as the data.
````

## Weighted RMSD
Instead of the standard RMSD, a weighted RMSD (wRMSD) may be used to prioritize 
certain datapoints over the others. Given a vector $\textbf{w}$ of $N$ weights $w_n$, one for each pair of 
points $\hat{\textbf{p}}_n$ and $\textbf{p}_n$ in $\hat{\textbf{P}}$ and $\textbf{P}$, respectively, 
the wRMSD is defined as {cite}`RapidRMSD`:

```{math}
:label: wRMSD
\begin{align}
\textrm{wRMSD}(\hat{\textbf{P}}, \textbf{P}, \textbf{w})
&= \sqrt{\frac{1}{\sum_{n = 1}^{N}w_n}\sum_{n = 1}^{N}w_n\|\hat{\textbf{p}}_n - \textbf{p}_n\|^2} \\
&= \sqrt{\frac{1}{\sum_{n = 1}^{N}w_n}\sum_{n = 1}^{N}w_n\sum_{k = 1}^{K}\left(\hat{p}_{nk} - p_{nk}\right)^2}
\end{align}
```

````{note}
Slightly different formulations of wRMSD are sometimes encountered in the literature. For example,
in a publication on using a Gaussian-weighted RMSD calculation for alignment of flexible proteins,{cite}`GaussianWRMSD`
the authors found that separating the weighted sum and the sum of all weightes was more suitable for their
application; that is, the wRMSD was defined as:

```{math}
\textrm{wRMSD}(\hat{\textbf{P}}, \textbf{P}, \textbf{w})
\sqrt{\frac{1}{N}\sum_{n = 1}^{N}w_n\|\hat{\textbf{p}}_n - \textbf{p}_n\|^2}
```
As a result, the absolute values of weights have an effect on the wRMSD value in this definition; that is,
scaling all weights by a constant factor also scales the wRMSD value by the square root of that factor.
Notice that in most alternative definitions, the values only differ from those obtained
by the eq. {eq}`wRMSD` by a constant factor. Therefore, as long as calculations are consistently done
with one definition, and it is clear how to interpret the values based on how the wRMSD was defined,
any one definition can be used. Nevertheless, most implementations, including the implementation in openCADD, 
follow the definition in eq. {eq}`wRMSD`
(cf. wRMSD implemented in [MDAnalysis](https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsd.html#Background)
[[API](https://docs.mdanalysis.org/stable/documentation_pages/analysis/rms.html#MDAnalysis.analysis.rms.rmsd),
[source code](https://github.com/MDAnalysis/mdanalysis/blob/develop/package/MDAnalysis/analysis/rms.py)]). 
````


```{hint}
In calculating the RMSD between two conformations of a molecule, the atomic mass
of each atom can be used as weight. Most molecules are mainly made of hydrogen, carbon, nitrogen, and oxygen
atoms. Carbon, nitrogen and oxygen have similar masses, whereas the mass of hydrogen is more than one order of
magnitude smaller than the others. Therefore, using atomic masses as weights for calculating the wRMSD
has the main effect that deviations in positions of hydrogen atoms are about 10 times less severly
penalized compared to other atoms.  
```

This definition can also be generalized to allow for different weights for different dimensions of the
points in the dataset. This may not be sensible for the case of molecular conformations, 
where the points are the coordinates of atoms, since there is in general no difference between
distances along x-, y- and z-axes. However, when the points describe different features
in an abstract feature spaces, then it may be useful to be able to weight each feature separately.
More generally, for each individual dimension of each individual point in the dataset, a different weight
may be assigned in a matrix of weights $\textbf{W}$, with the same shape as the dataset ($N\times K$):

```{math}
:label: W
\textbf{W}_{N\times K} 
= \begin{bmatrix}
    \hat{\textbf{w}}_1 \\
    \vdots \\
    \hat{\textbf{w}}_{N}
    \end{bmatrix}
= \begin{bmatrix}
    \hat{w}_{11} & \dots & \hat{w}_{1K} \\
    \vdots & \ddots & \vdots \\
    \hat{w}_{N1} & \dots & \hat{w}_{NK}
    \end{bmatrix}
```

The wRMSD is then defined as (notice that the sum of all weights is 
divided by the dimension of the points, in order to remain consistent 
with the prior definitions):
```{math}
\begin{align}
\textrm{wRMSD}(\hat{\textbf{P}}, \textbf{P}, \textbf{W})
&= \sqrt{\frac{1}{\frac{\sum_{n = 1}^{N}\sum_{k = 1}^{K}w_{nk}}{K}}\sum_{n = 1}^{N}\sum_{k = 1}^{K}w_{nk}\left(\hat{p}_{nk} - p_{nk}\right)^2} \\ 
&= \sqrt{\frac{K}{\sum_{n = 1}^{N}\sum_{k = 1}^{K}w_{nk}}\sum_{n = 1}^{N}\sum_{k = 1}^{K}w_{nk}\left(\hat{p}_{nk} - p_{nk}\right)^2} 
\end{align}
```



