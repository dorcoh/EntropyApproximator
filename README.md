# Info

<p>A supplement for the book <a href="https://mitpress.mit.edu/books/perturbations-optimization-and-statistics">&quot;Perturbations, Optimization, and Statistics&quot;</a> (edited by Tamir Hazan, George Papandreou and Daniel Tarlow) - Chapter 9: &quot;Probabilistic Inference by Hashing and Optimization&quot; by Stefano Ermon</p>

A talk about Chapter 9 is also available at `Seminar.pdf`

<p><b>Notice:</b> The tools attached here (STS, WISH) are under copyright of Stefano Ermon , his README's are inside the directories</p>

<p>For the latest versions you can check his page: http://www.cs.cornell.edu/~ermonste</p>

# Patch

<p>In my <a href="https://dorcoh.github.io/entropy-patch/">blog post</a> I describe a patch to Ermon's SampleTreeSearch tool, to compute the approximated entropy (the average degree of freedom for it's variables) of a CNF formula.</p>

<p>Applying the patch could be done in the following way (tested under Linux Ubuntu):</p> 

### Patch main (can skip this and instead download and compile STS-Entropy):

```
[user@localhost ~]$ patch Main.cc -i entropy.patch -o newMain.cc
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
Average Entropy: 0.415007
```

# References

### STS
Stefano Ermon, Carla Gomes, and Bart Selman.

Uniform Solution Sampling Using a Constraint Solver As an Oracle.

UAI-12. In Proc. 28th Conference on Uncertainty in Artificial Intelligence, August 2012.

### 
