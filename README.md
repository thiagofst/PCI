# Propagation of Confidence Interval (PCI)

PCI is a package for handling measurement values with asymmetric uncertainties given by their confidence intervals.

Measured values and its respective confidence intervals are sampled with a pre-defined likelihood function, called "Variable Width Gaussian", 
proposed by R. Barlow's 2004 paper "Asymmetric Statistical Errors".


The operations are performed by Monte Carlo simulation, simple operations (+,-,*,/ ,**) are straightforward and the package has support for generic functions to be created by the user. For now the code supports 68%, 90% and 95% confidence intervals:

                                                  Confidence	Δχ2     ΔlogL
                                                      68%	       1.00     -0.50
                                                      90%	       2.71     -1.36
                                                      95%	       4.00     -2.00



## Installation

Use the package manager [pip](https://pip.pypa.io/en/stable/) to install PCI.

```bash
pip install git+https://github.com/muryelgp/PCI.git
```

## How to use the code

```bash
from pci import AssymetricMeasurment as asyme
#variable(measurment, err_n, err_p, confidence)
a=asymed(10, 1, 1.5, confidence=68)
b=asymed(30, 3, 3.5, confidence=68)
c = a + b 
```
