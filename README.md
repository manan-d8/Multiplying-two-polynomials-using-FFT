# Multiplying-two-polynomials-using-FFT
This is Assignment 1 of CS6013: Advanced Data Structures and Algorithms [IITH] - Multiplying two polynomials


-------------------------------------------------------------------
# Programming Assignment 1 
### Name        : Manan Darji
### Roll Number : ************
### Subject     : Advance Data Structure and Algorithm (CS6013)
### Topic       : Multiplying Two Polynomials
-------------------------------------------------------------------

### Flow Of Code:
* I Have Used **Fast Fourier Transform** to Implement problem of Multiplying Two Polynomials.
* Code will give you both output.
  1) For naive approach
  2) For FFT approach

<br/>

1) **For naive approach**
   * For Naive Approach we'll just multiply each term of first polynomial with each term of second polynomial.
   * This will take O(N^2) Time. 

<br/>

2) **For FFT approach**
   * Using Fast Fourier Transform will Give Us Time complexity as O(NLog(N)).
   * first we'll convert both polynomial to point value representation in O(NLog(N)) using Fast Fourier Transform.
   * then we'll  multiply both of them using simple multiplication of complex number in O(N).
   * After which we'll convert resulting point value representation to polynomial coefficient representation in O(NLog(N)) using Inverse Fast Fourier Transform.
   * At end we normalized the output of IFFT.





<div style="page-break-after: always"></div>

### How output of code looks like
```
----------------------------------------------------------------------------------
Please enter degree of Polynomial 1 : 4
Please enter 5 coefficients for polynomial 1 in the increasing order of the degree of the monomials
5 4 8 7 9
----------------------------------------------------------------------------------
Please enter degree of Polynomial 2 : 3
Please enter 4 coefficients for polynomial 2 in the increasing order of the degree of the monomials
2 5 4 8
----------------------------------------------------------------------------------
polynomial 1 : 9x*4 + 7x*3 + 8x*2 + 4x + 5
polynomial 2 : 8x*3 + 4x*2 + 5x + 2
----------------------------------------------------------------------------------
polynomial Product (Naive) : 72x*7 + 92x*6 + 137x*5 + 117x*4 + 110x*3 + 56x*2 + 33x + 10
----------------------------------------------------------------------------------
polynomial Product ( FFT ) : 72x*7 + 92x*6 + 137x*5 + 117x*4 + 110x*3 + 56x*2 + 33x + 10
----------------------------------------------------------------------------------
Both Naive And FFT Output Are Same, All Good!
----------------------------------------------------------------------------------
```
