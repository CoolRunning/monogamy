## Quantum Monogamy ##

### About ###

This folder contains a small python-module that allows you to discover multi-partite quantum-correlations. It was created as a research project for the lecture "Quantum Computation and Communication" given by Stephanie Wehner 2015 at TU Delft.

### Introduction ###

Quantum Monogamy is a concept from Quantum Information theory. It describes the relation of different parties that may share (entangled) quantum systems. In particular, it can be shown that if two parties, Alice and Bob, share a maximally entangled state (e.g.: EPR-pair), there can be no third party (Charlie) be entangled with either Alice's or Bob's subsystem. Alice and Bob have stablished a monogamous correlation with each other, that limits the access on information for Charlie. However, if the correlation between Alice and Bob is not perfect but subject to noise, Charlie might share state information with Alice and could gain information about measurement outcomes between Alice and Bob.

In order to design robust quantum cryptographic protocols it is of high importance to investigate to what extent quantum states can be correlated between multiple parties. 

### Installing the module ###

The module is written in [Python3](https://www.python.org). It depends on [NumPy](http://www.numpy.org/) and [CVXPY](http://cvxpy.readthedocs.org/en/latest/index.html). If you want to use it, you first have to set up a working python3 environment with CVXPY on your system. Installation instructions are found [here](http://cvxpy.readthedocs.org/en/latest/install/index.html).

Once your python-environment runs CVXPY, clone the github-repository. If you start a python-shell in the cloned folder (recommended: [IPython](http://ipython.org/)) you will be able to import and access the module for example by

```python
import monogamy as mg
```
 
The interface of monogamy is documented in the examples below.

### Using the module ###

The following IPython-notebook illustrate how to use the module. Once you downloaded the IPython-notebooks, you can experiment with them locally and change some values to see what happens.

* [Construction of global states](http://nbviewer.ipython.org/github/CoolRunning/monogamy/blob/master/global.ipynb)

* [Usage of the API](http://nbviewer.ipython.org/github/CoolRunning/monogamy/blob/master/mono.ipynb)

### Monogamy of Werner States ###

As a research project, the module was supposed to be used to analyse a family of quantum states called [Werner States](http://en.wikipedia.org/wiki/Werner_state). The following IPython-notebook summarizes the task and the obtained results, using the monogamy module.

* [Werner State monogamy](http://nbviewer.ipython.org/github/CoolRunning/monogamy/blob/master/werner_monogamy.ipynb)

### License ###

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
