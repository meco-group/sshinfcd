# sshinfcd - State-space $\mathcal{H}_\infty$ control design for LTI systems

`sshinfcd` is a MATLAB package facilitating multi-objective $\mathcal{H}_\infty$ controller synthesis with $D$-stability constraints. It is a relatively straightforward implementation of a combination of the LMI formulations discussed in

* H. Köroğlu. ['$\mathcal{H}_\infty$ synthesis with unstable weighting filters: An LMI solution.'](http://dx.doi.org/10.1109/CDC.2013.6760244) 52nd IEEE Conference on Decision and Control, pp. 2429-2434, 2013.
* M. Chilali and P. Gahinet. ['$\mathcal{H}_\infty$ design with pole placement constraints: an LMI approach.'](http://dx.doi.org/10.1109/9.486637) IEEE Transactions on Automatic Control, vol. 41, no. 3, pp. 358-367, 1996. 

The work presented in the first reference has also been extended for descriptor realizations and impulsive weighting filters (see, e.g., [Feng & Yagoubi, 2016](https://doi.org/10.1016/j.automatica.2016.01.028) or [Feng et al., 2020](https://doi.org/10.1080/00207179.2018.1562223)). A general implementation of those techniques that is numerically sufficiently reliable for practical problems, however, turns out to be challenging. In addition, it has been shown that, if weighting filters are not impulsive, the problem can always be transformed into a format that can be solved with `sshinfcd` using regularization ([Bunjaku et al. 2018](https://doi.org/10.1109/ICCA.2018.8444268), [Wang et al., 1987](https://doi.org/10.1016/0167-6911(87)90028-4)). 

It is important to note that the The Mathworks' Control Systems Toolbox for use with MATLAB also allows the same problem type to be solved using `systune`. In no way, this package aims at *competing* with these existing tools. Instead, we consider it to be *complementary*. Whereas `systune` is a lot more flexible in terms of specifications and structural constraints, it uses a nonconvex reformulation of the design problem, implying that the solution may heavily depend on the initial guess. When many initial guesses are required, the procedure also tends to become rather slow. In contrast, `sshinfcd` uses a classic convex approach by reformulating the problem as a semi-definite program (SDP) through linear matrix inequalities (LMIs). As long as the order of the system is moderate, modern SDP solvers can quickly return accurate results. The associated downside of an SDP reformulation is the conservatism that is introduced to render the formulation convex. Both solution methods are useful for practicing control engineers, depending on the needs of their specific problem at hand.

## Prerequisites

blabla

## Usage

blabla

### Basic functionality

blabla

### Options

blabla

## A note on LCToolbox

blabla



