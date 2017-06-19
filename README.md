<p>A supplement for the book <a href="https://mitpress.mit.edu/books/perturbations-optimization-and-statistics">&quot;Perturbations, Optimization, and Statistics&quot;</a> (edited by Tamir Hazan, George Papandreou and Daniel Tarlow) - Chapter 9: &quot;Probabilistic Inference by Hashing and Optimization&quot; by Stefano Ermon</p>

<p><b>Notice:</b> The tools described here are under copyright of Stefano Ermon , his README's are inside the directories</p>

<p>For the latest versions you can check his page: http://www.cs.cornell.edu/~ermonste</p>

<p>In my <a href="https://dorcoh.github.io/entropy-patch/">blog post</a> I describe a patch to Ermon's SampleTreeSearch tool, to compute the approximated entropy (the average degree of freedom for it's variables) of a CNF formula.</p>

<p>Applying the patch could be done in the following way (under Linux):</p>

```
[user@localhost ~]$ patch originalfile -i patchfile.patch -o updatedfile
```

A presentation of mine is also available `Seminar.pdf`