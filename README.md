# Message passing for spectral density of sparse matrices

This is an implementation of the message passing algorithm to compute the spectral density of sparse matrices, introduced in "*[Message passing on networks with loops](https://arxiv.org/abs/1907.08252)*".


## Compiling and running

Run `make` to compile.

Options:<br>
+  -i &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; input file<br>
+  -o &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; output file<br>
+  -r &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; neighborhood level, r=0,1,2<br>
+  -z &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; points to evaluate rho(z), format:  x_min,x_max,number_of_points,Im(z)<br>
 
For example,
```
./spectrum -i example/laplacian.txt -o out.txt -r 2 -z 0,10,50,0.1
```
will run the algorithm on the laplacian matrix saved in laplacian.txt, and output to out.txt.  It will run using the r=2 neighborhood, and will compute the value of rho(x) at 50 points between 0 and 10. The broadening parameter, eta, will be set to 0.1.

Tested using gcc 9.2.0.

## Reference

If you use this code, please consider citing:

"[*Message passing on networks with loops*](https://arxiv.org/abs/1907.08252)"<br/>
George T. Cantwell and M. E. J. Newman<br/>
Proc. Natl. Acad. Sci. USA (in press) <br/>
