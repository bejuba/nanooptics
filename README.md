# nanooptics
Python module for data evaluation and presentation developed in the nanooptics group at the Humboldt University Berlin.

## installation
The correlator makes heavy use of cython. Therefore those modules have to be compiled.
The setup.py should take care of this when installing the module with pip or by executing the setup.py.

Therefore simply cloning and installing should just work


```
git clone https://github.com/bejuba/nanooptics.git

cd nanooptics

python setup.py install
```

or using pip


```
pip install https://github.com/bejuba/nanooptics/archive/master.zip

```

## fixing cython compiler errors 

if your systems gcc is too new (e.g current arch linux installations), it might not work with the ld required for cython.
if you are using a anaconda installation you can ensure that the anaconda compilers are getting used instead by installing them with:

```
conda install gxx_linux-64
```
see
https://github.com/Anaconda-Platform/anaconda-project/issues/183

## Examples

the examples folder contains jupyter-notebooks and a tcspc smalll data set in the form of a .pt2 file measured with a picoharp as a starting point for common data evaluation schemes.
