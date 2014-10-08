Analytical Inverse Optimization (ANIO) Toolbox for MATLAB
=========================================================

Distributes under GNU General Public Licence, see `LICENSE` file.

The toolbox implements the basic scripts for ANIO analysis. In the
current version it offers two algorithms for estimating the cost
function: `anio` and `anioangle`. These algorithms are described in
details in [Terekhov et al 2010] [1] and [Terekhov & Zatsiorsky 2011]
[2]. These papers also describe the conditions on the data under which
the algorithms can be used at all.

Preamble
--------

ANIO aims at finding the conditions under which a cost function can be
estimated unambiguously from experimental data. It, thus, assumes that
there is a certain data set `X`, which is obtained as a result of
minimization of a certain cost function

        J(x) --> min
	
under constraints

        C x = b
	
It is assumed that there is a data set `X` of vectors `x`, obtained
for different task conditions `b` and the same matrix of task
constraints `C`.

The theorem of uniqueness states that the cost function `J(x)` can be
determined from the data set `X` if:

1.  The cost function J(x) is linear additive, meaning that

    `J(x1,...,xn) = g1(x1) + ... + gn(xn)`

    where g1, ..., gn are unknown functions of a scalar argument and
    x1, ..., xn are the components of the vector x.

2.  The data is corrected for a certain volume in the space of b. For
    example if b is 2-dimensional than the data must be collected for
    a certain area in this 2d space, not for a line. Note, that the
    dimensionallity of b needs to be **at least two**.
    
3.  The constraints are 'non-splittable', meaning that the matrix
    `C'*inv(C*C')*C cannot be made block-diagonal by simultaneous
    swapping of columns and rows. For example, a matrix

    ` 1 0 2 `
    
    ` 0 3 0 `
    
    ` 2 0 4 `

    is splittable because simultaneous swapping the first and the
    second rows and columns will transform it into

    ` 3 0 0 `
    
    ` 0 1 2 `
    
    ` 0 2 4 `

    which is block-diagonal.

The cost function `J(x)` can be determined 'almost unambiguously',
which means that if the data set was produced using a function `J(x)`,
the estimated function `Ĵ(x)` is equal to

        Ĵ(x) = r*J(x) + q'*C*x + const
	
where `r` is a positive scalar, `q` is a vector of the same dimension
as `b` and `const` is a scalar.

Demonstration problem
---------------------

As a test example consider a task in which an imaginary subject is
instructed to produce a certain force and moment of force by pressing
with four fingers on force sensors. Spacing between the sensors is 3
cm.

In the experiment the force changes from 10 to 15 N with the step of 1
N and the moment of force changes from -0.1 Nm to 0.1 Nm, with a step
of 0.05 Nm. The imaginary subject performs the task for every
combination of force and moment of force. It is programmed to choose
forces that minimize the cost function

        J(x1,x2,x3,x4) = 1/2 * (k1 x1^2 + k2 x2^2 + k3 x3^2 + k4 x4^2) +
	                      + v1 x1   + v2 x2   + v3 x3   + v4 x4

subject to the constraints

        C x = b

The data for the imaginary subject is saved in `example1.mat`. The
file contains:

-   the matrix of constraints `C`, in which the first row corresponds to
    the constraints on the total force and the seconds corresponds to
    the constraints on the total moment of force;

-   the matrix of the target values `b`, the first column corresponds to
    the target force, the second column corresponds to the target
    moment of force; each row corresponds to an experimental condition;

-   the matrix of 'clean' data `X0`, in which each row
    corresponds to the same experimental condition as in `b`;

-   the matrix of the 'measured' data `X` corrupted by the motor
    variability and the measurement noise, which were simulated as a
    Gaussian noise with the standard deviation of 0.3 N;

-   the 'true' coefficients `w0`, where the first column corresponds to
    the linear coefficients and the second column corresponds to
    the quadratic coefficients.

Example of the scripts usage
----------------------------

Change your Matlab current folder to the one where you put ANIO
toolbox files and load the file `example1.mat` in your Matlab workshop
by typing in the command line

        load('example1.mat');

You can take a look at the data by plotting its 3d projections

        plot(X(:,1), X(:,2), X(:,3), '.')

to see if it tends to form a 2d surface.

This can be further reassured by running PCA and checking the
percentage of accounted variance:

        L = eig(cov(X)); L'/sum(L)*100

which should give

        78.0733    18.8187    1.8441    1.2639

So that the two main PCs account or about 97% of the total
variance. Moreover, the percentage of variance accounted for by the
second major PC is dramatically higher than that by the third PC,
confirming that the data are distributed along a plane.

The fact that the data are distributed in a plane suggests that the
cost function must be the second order polynomial of the forces. For
details see [Terekhov et al 2010] [1]. In this case `anioangle` can be used
to determine the cost function coefficients. Type

        [w, ang] = anioangle(C, X);

to get the cost function coefficients `w` and the dihedral angle `ang`
between the plane defined by the first 2 PCs and the optimal plane.

        ang * 180/pi

should return 2.9043. Small value of the dihedral angle, together with
good PCA scores, indicates that the data can indeed be fitted by a
cost function additive with respect to the finger forces.

The estimated coefficients are returned in `w`, where the first column
contains first order coefficients and the second column contains
second order coefficients.

Typing just

        w

should produce

         0.0208    0.5949

         0.0952    0.2844

        -0.2530    0.4160

         0.1369    0.6262

The first order coefficients are significantly different from the true
ones stored in `w0`, but this is precisely what is predicted by the
theorem of uniqueness in [Terekhov et al 2010] [1]. The second order
coefficients are much closer to their true values.

The cost function coefficients can also be estimated using a more
general method, described in [Terekhov & Zatsiorsky] [2]. To do that,
run

        w = anio(X, C);

The return values should be

         0.0030    0.6331

         0.0093    0.3000

        -0.0274    0.3603

         0.0152    0.6160

Again, the second order coefficients are rather close to their true values.

For more information contact Alexander Terekhov (avterekhov AT
gmail.com) or visit the Facebook page for ANIO:
https://www.facebook.com/anio.method

[1]: http://link.springer.com/article/10.1007/s00285-009-0306-3
[2]: http://link.springer.com/article/10.1007/s00422-011-0421-2
