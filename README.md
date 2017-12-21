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
Output file: t2307.cnf.entropy.out
```
Output file: <a href="https://github.com/dorcoh/entropyApproximator/blob/master/STS-Entropy/core/t2307.cnf.entropy.out">t2307.cnf.entropy.out</a>

### Args:
As described in `STS --help`:
```
  -nsamples     = <int32>  [   0 .. 300000000] (default: 10)
  -k            = <int32>  [   0 .. 100000000] (default: 50)
```
nsamples - Number of sampling iterations
k - Number of samples per level (the higher the value, the more uniform the solutions are)

# Sanity check

Let's examine a simple formula with 3 variables:

```
p cnf 3 3
1  2  3  0
1  -2  3  0
1  2  -3  0
```

We can write it down and see that there are total 5 solutions:
` { (1,0,0), (1,0,1), (1,1,0), (1,1,1), (0,1,1) } `

Ratios of each literal (number of times it appears in solutions):
` r(1) = 4/5, r(-1) = 1/5, r(2) = 3/5, r(-2) = 2/5, r(3) = 3/5, r(-3) = 2/5 `

In this scenario of tiny formula the approximator sampled 50 solutions. Some of the solutions are identical of course, but usually it isn't the case where we handle larger formulas. In particular the samples should be sampled uniformly (randomly). Let's take a look of it's output:

```
Var,TotalSols,PosLitSols,NegLitSols,EntropyShan
1,50,0.800000,0.200000,0.721928
2,50,0.600000,0.400000,0.970951
3,50,0.600000,0.400000,0.970951
#Estimated entropy: 0.887943
```

The ratios (PosLitSols/NegLitSols) converged <b>exactly</b> to the right values.

# References

### STS
Stefano Ermon, Carla Gomes, and Bart Selman.

Uniform Solution Sampling Using a Constraint Solver As an Oracle.

UAI-12. In Proc. 28th Conference on Uncertainty in Artificial Intelligence, August 2012.

### About entropy property of CNF formulas
Dor Cohen, Ofer Strichman

The impact of Entropy and Solution Density on selected SAT heuristics

abs/1706.05637, arXiv pre-print, June 2017