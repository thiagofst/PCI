# Propagation of Confidence Intervals (PCI)

PCI is a package for handling measured values with asymmetric uncertainties given by their confidence level in which the probability distribution is unknown.

The values and their respective confidence intervals are sampled with a pre-defined likelihood function, called "Variable Width Gaussian", 
proposed in [R. Barlow's 2004 paper "Asymmetric Statistical Errors"](https://arxiv.org/abs/physics/0406120).


The operations are performed by Monte Carlo simulations, simple operations (+, -, *, /, **) are straightforward and the package has support for generic functions to be created by the user. For now, the code supports 68%, 90%, and 95% confidence levels:

                                              Confidence	Δχ2     ΔlogL
                                                  68%     1.00     -0.50
                                                  90%     2.71     -1.36
                                                  95%     4.00     -2.00



## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install PCI.

```bash
pip install git+https://github.com/muryelgp/PCI.git
```

## How to use the code

```bash
from pci import AssymetricMeasurment as asyme
```
A measurment with a nominal value of 10 and errors +1 and -0.5 at 68% confidence level or: 

![equation](https://latex.codecogs.com/gif.latex?a&space;=&space;10^{&plus;1.0}_{-0.5})

can be instanciated by:

```bash
a=asymed(10, 0.5, 1, confidence=68)
```
in the same way:

![equation](https://latex.codecogs.com/gif.latex?b&space;=&space;30^{&plus;3.0}_{-3.5})


```bash
b=asymed(30, 3.5, 3, confidence=68)
```
Simple operations are straightforward, for example: 

```bash
c = a + b 
```

Check out this [jupyter notebook](https://github.com/muryelgp/PCI/blob/master/pci/How_to.ipynb) for details and more interesting examples.
