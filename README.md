# The Error Function
I explored several different algorithms for the computation of the error
function. The version included in `libtmpl` had the best speed while attaining
the appropriate precision. The rejected algorithms are collected here for
future study.

Note that while some algorithms are much faster than the one in `libtmpl`, they
only achieve a few decimals of accuracy, and hence shouldn't be used in a
`libm` implementation.

## Building
The only dependency is [`libtmpl`](https://github.com/ryanmaguire/libtmpl).
Build and install this first.
You can then compile for `C` use using the Makefile.
```sh
make
```
You can choose compiler with the `CC` option.
```sh
make CC=clang
```
And add extra flags if you want.
```sh
make CC=clang Extra_Flags=-Weverythin
```

A small Python wrapper is available with `numpy` support so you can plot the
functions with `matplotlib`. To build this make sure `numpy` and `python3`
are installed. After building and installing `libtmpl, type
```sh
python3 setup.py build_ext --inplace
```
