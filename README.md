# Message passing for spectral density of sparse matrices

This is an implementation of the message passing algorithm to compute the spectral density of sparse matrices, introduced in *[Message passing on networks with loops](https://doi.org/10.1073/pnas.1914893116)*.


## Compiling and running

Run `make` to compile.

Options are passed on command line, in the following order:
./spectrum input_matrix.txt r eta x_min x_max num_points
where:
input_graph.txt is the sparse input matrix
r is the neighborhood level (r=0,1,2,3...)
eta is the broadening parameter
spectral density evaluated between x_min and x_max and num_points number of points.

For example,
```
./spectrum example/laplacian.txt  2 0.1 0 10 50
```
will run the algorithm on the laplacian matrix saved in laplacian.txt.  It will run using the r=2 neighborhood, and will compute the value of rho(x) at 50 points between 0 and 10. The broadening parameter, eta, will be set to 0.1.

By default this software uses OpenMP to run in parallel.  If this isn't available it can be compiled to run on a single thread with `make single_thread`.


## Reference

If you use this algorithm please cite:

[*Message passing on networks with loops*](https://doi.org/10.1073/pnas.1914893116)<br/>
George T. Cantwell and M. E. J. Newman<br/>
Proc. Natl. Acad. Sci. USA 116 (2019)<br/>
