Approximate a CNF formula solutions entropy (defined in paper referenced below), by using hashing & optimization technique as implemented in STS SAT model counter (reference below).

# Info

<p>A supplement for the book <a href="https://mitpress.mit.edu/books/perturbations-optimization-and-statistics">&quot;Perturbations, Optimization, and Statistics&quot;</a> (edited by Tamir Hazan, George Papandreou and Daniel Tarlow) - Chapter 9: &quot;Probabilistic Inference by Hashing and Optimization&quot; by Stefano Ermon</p>

A talk about Chapter 9 is also available at `Seminar.pdf`

<p><b>Notice:</b> The tools attached here (STS, WISH) are under copyright of Stefano Ermon , his README's are inside the directories</p>

<p>For the latest versions you can check his page: http://www.cs.cornell.edu/~ermonste</p>

# Patch

<p>In my <a href="https://dorcoh.github.io/entropy-patch/">blog post</a> I describe a patch to Ermon's SampleTreeSearch tool, to compute the approximated entropy (the average degree of freedom for it's variables) of a CNF formula.</p>

<p>Applying the patch could be done in the following way (tested under Linux Ubuntu):</p> 

### Patch main

(can skip this and instead download and compile STS-Entropy)

```
[user@localhost ~]$ patch Main.cc -i entropy.patch -o Main.cc
```

### Compile:
```
cd <STS-DIR>
export MROOT=$PWD
cd core
make

```

### Test:
```
./STS ../test/t2307.cnf
```

### Output:
```
Different : 457
Chi-square : 39.260000
Estimated log-z: 7.666618
Estimated Z: 2.135845e+03
Estimated entropy: 0.472847
Output file: <a href="https://github.com/dorcoh/entropyApproximator/blob/master/STS-Entropy/core/t2307.cnf.entropy.out">t2307.cnf.entropy.out</a>
```

## Sanity check

Let's examine a simple formula with 3 variables:

$` x_{1} , x_{2}, x_{3} `$

Constraints (CNF formula):

$ \varphi = (x_{1} \lor x_{2} \lor x_{3} ) $

# References

### STS
Stefano Ermon, Carla Gomes, and Bart Selman.

Uniform Solution Sampling Using a Constraint Solver As an Oracle.

UAI-12. In Proc. 28th Conference on Uncertainty in Artificial Intelligence, August 2012.

### About entropy property of CNF formulas
Dor Cohen, Ofer Strichman

The impact of Entropy and Solution Density on selected SAT heuristics

abs/1706.05637, arXiv pre-print, June 2017