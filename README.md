## Quantum Monogamy ##

### About ###

This folder contains a small python-module that allows you to discover multi-partite quantum-correlations. It was created as a research project for the lecture "Quantum Computation and Communication" given by Stephanie Wehner 2015.

### Introduction ###

Quantum Monogamy is a concept from Quantum information theory. It describes the relation of different parties that may share (entangled) quantum systems. In particular it can be shown that if two parties, Alice and Bob, share a maximally entangled state (EPR-pair), there can be no third party (Charlie) be entangled with either Alice's or Bob's subsystem. Alice and Bob have stablished a monogamous correlation with each other, that limits the access on information from Charlie. However, if the correlation between Alice and Bob is not perfect but subject to noise, Charlie might share state information with Alice and could gain information about measurement outcomes between Alice and Bob.

In order to design robust quantum cryptography protocols it is thus of importance to investigate to which extend quantum states can be correlated between multiple parties. 

### Installing the Module ###

The module is written in requires [Python3](https://www.python.org). It depends on [NumPy](http://www.numpy.org/) and [CVXPY](http://cvxpy.readthedocs.org/en/latest/index.html). If you want to use it, you first have to set up a working python3 environment with CVXPY on your system. Installation instruction are found [here](http://cvxpy.readthedocs.org/en/latest/install/index.html).

Once your python-environment runs CVXPY, create an empty folder and download monogamy.py from this repository. If you start a python-shell in this folder (recommended: [IPython](http://ipython.org/)) you will be able to import and access the module simply by

```python
import monogamy as mg
```
 
Follow some of the examples to get familiar with the interface.

### Using the module ###

### Monogamy of Werner States ###

### License ###

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
