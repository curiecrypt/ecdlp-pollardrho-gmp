# ecdlp-pollardrho-gmp
This C code is the GMP implementation of Pollard's Rho Attack on Elliptic Curve Discrete Logarithm Problem (ECDLP).

## Getting Started

In order to compile this code, a C compiler and [GMP Library](https://gmplib.org/)  (The GNU Multiple Precision Arithmetic Library) is required. gcc version 6.3.0 and GMP 6.1.2. have been used for testing.

[Magma](http://magma.maths.usyd.edu.au/magma/) Computational Algebra System is used as an intermediate step between theory and programming. Also, the inputs of C code are obtained from Magma (See Testing section). You can purchase the software or can use the [Online Calculator](http://magma.maths.usyd.edu.au/calc/). Open-source lovers can use [SageMath](http://www.sagemath.org/) as well but this repository only contains Magma scripts.


## Testing

First step is to create elliptic curve and choose the random variables in Magma. The following functions will be needed:

```
sum := function(P,Q,A,E)

	x1 := P[1];
	x2 := Q[1];
	y1 := P[2];
	y2 := Q[2];

	if P eq E!0 then return Q;
	elif Q eq E!0 then return P;
	elif x1 eq x2 then
		if y1 ne y2 then return E!0;
		elif y2 eq 0 then return E!0;
		else
			x3 := ((3*x1^2 + A )/(2*y1))^2 - 2*x1 ;
			y3 := (3*x1^2 + A )/(2*y1)*(x1-x3) - y1 ;
			return E![x3,y3];
		end if ;
	else
		x3 := ((y1- y2)/(x1-x2))^2 - x1- x2 ;
		y3 := (y1- y2)/(x1-x2) * (x1-x3) -y1 ;
		return E![x3,y3];
	end if ;

end function;

mult := function(k,P,E,A)
	M :=E!0;
	K := IntegerToSequence(k,2);

	for i in [1..#K] do
		if K[i] eq 1 then
			M := sum(M,P,A,E);
		end if;
		P := sum(P,P,A,E);
	end for;
	return M;
end function;
```

Then, you can obtain the inputs for attack:

```
p := NextPrime(Random([0..2^64]));
F := FiniteField(p);

repeat
	repeat
		A := Random(F);
		B := Random(F);
	until (4*A^3 +27*B^2) ne 0 ;
	E := EllipticCurve([A,B]);
until IsPrime(#E);

P := Random(E);
n := Order(P);
k := Random([0..n]);
Q := mult(k,P,E,A);

P; Q; n; k;
```
The full Magma code can be used for the verification of results of the C code.

Note: Magma Online Calculator has time limitation, thus Magma code cannot be run on Online calc. for 64 bit parameters.
