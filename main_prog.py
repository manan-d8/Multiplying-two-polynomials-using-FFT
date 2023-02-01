"""
-------------------------------------------------------------------
Programming Assignment 1 
-------------------------------------------------------------------
Name        : Manan Darji
Roll Number : CS22MTECH14004
Subject     : Advance Data Structure and Algorithm (CS6013)
Topic       : Multiplying Two Polynomials
-------------------------------------------------------------------
"""

# Imported math to use Cos and Sin Function For finding the roots.
import math


# Defined 2 Polynomials to take input | max size of degree = 15.
poly_1 = [0] * 16
poly_2 = [0] * 16

# Array to store naive multiplication.
naive_prod = [0] * 32

# Array to store FFT multiplication.
fft_prod = [0] * 32

# Array to store N nth roots of unity globally.
Roots = []


def naive_polynomial_multiplication(poly_1, poly_2):
    """This Function is Used for Naive Multiplication here we do Multiplication
    in O(n^2) Time, Just by multiplying every term of first polynomial with
    every term of second polynomial.
        It Will store Output in Global Array naive_prod.

    Args:
        poly_1 (int[]): First Polynomial
        poly_2 (int[]): Second Polynomial
    """
    for i, val1 in enumerate(poly_1):
        for j, val2 in enumerate(poly_2):
            if val1 * val2:
                naive_prod[i + j] += val1 * val2


def print_poly(Coefficient_Array):
    """This Function is used for printing the polynomial if natural format.
        ex: 1) 12x*2 + 17x + 6+
            2) 3x + 2
            3) - 55x*4 - 8x*3 + 7x*2 + 45x + 2

    Some feature :
        1) it will not print Coefficient if it's 1
            1x*3 --> x*3
        2) it will not print degree if it's 1
            5x*1 --> 5x
        3) it will only print sign of degree n term is it's negative
        4) it will not print terms with coefficient as 0

    Args:
        Coefficient_Array (int[]): Coefficient/vector representation of polynomial
    """
    First_Flag = False
    N = len(Coefficient_Array) - 1
    for i, val in enumerate(reversed(Coefficient_Array)):
        if First_Flag and val > 0:
            print(" + ", end="")
        elif val < 0:
            print(" - ", end="")
        if val != 0:
            if not N - i:
                print(abs(val), end="")
            else:
                print(
                    (abs(val), "")[abs(val) < 2],
                    "x",
                    ("*" + str(N - i), "")[(N - i) < 2],
                    sep="",
                    end="",
                )
            First_Flag = True
    print()


def add_complex_num(a, ai, b, bi):
    """This Function is used for summation of 2 complex number.

    Args:
        a  (Double): Real part of first complex number
        ai (Double): Imaginary part of first complex number
        b  (Double): Real part of second complex number
        bi (Double): Imaginary part of second complex number

    Returns:
        (Tuple): This function returns tuple of size 2 which is like (real,imag)
    """
    sum = (a + b, ai + bi)
    return sum


def multiply_complex_num(a, ai, b, bi):
    """This Function is used for Multiplication of 2 complex number.

    Args:
        a  (Double): Real part of first complex number
        ai (Double): Imaginary part of first complex number
        b  (Double): Real part of second complex number
        bi (Double): Imaginary part of second complex number

    Returns:
        (Tuple): This function returns tuple of size 2 which is like (real,imag)
    """
    mult = (((a * b) - (ai * bi)), ((a * bi) + (ai * b)))
    return mult


def find_N(deg_1, deg_2):
    """This Function Will return next number which is grater than
    deg1 + deg2 + 1 and it's in power of 2.

    Args:
        deg_1 (int): Degree of first  polynomial
        deg_2 (int): Degree of second polynomial

    Returns:
        int : returns number in power of 2
    """
    N, B = deg_1 + deg_2 + 1, 1
    if N and not (N & (N - 1)):
        return N
    while B < N:
        B <<= 1
    return B


def find_complex_roots(N):
    """This function will generate N nth roots of unity and will
    store it in the global array called roots.
        Here will always generate roots in power of 2 so N is in
    power of 2.

    Args:
        N (int): Number of roots to generate.
    """
    for i in range(N):
        Roots.append(
            ((math.cos((2 * math.pi * i) / N)), (math.sin((2 * math.pi * i) / N)))
        )


def Eval(POLY, IFFT):
    """This is the Recursive function which we use to do FFT and IFFT, basically based on
    the value of IFFT flag we implement FFT or IFFT.
        So we use this function to convert our vector representation into point value 
    representation and vice versa. for coefficient to point value representation conversion 
    we use FFT and for Point value to Coefficient we use IFFT.
        This Function will work in O(nLog(n)) Time.
        
    Args:
        POLY (array(tuple)):  for FFT it's array of tuple[2] where each tuple[2] represent complex number.
        IFFT (Boolean): Flag which represents that is it IFFT or FTT.

    Returns:
        array of tuple[2] where each tuple[2] represent complex number.
    """

    # This is length of input array
    N = len(POLY)

    # Base Condition for Recursion
    if N == 1:
        if IFFT:
            return POLY
        return [(POLY[0], 0)]
        # It return the coefficient in complex number / tuple format.

    # Here we Divide our array into two part in order to Implement Recursive Algorithm.
    POLYEven = []
    POLYOdd = []

    for i, val in enumerate(POLY):
        if i % 2 == 0:
            POLYEven.append(val)
        else:
            POLYOdd.append(val)

    # This is the recursive call for both new polynomials of size n/2 with even and odd terms.
    Eeven, Eodd = Eval(POLYEven, IFFT), Eval(POLYOdd, IFFT)

    # In this array we merge two recursive call output.
    EvalReturn = [(0, 0)] * N

    # So this loop will run for N/2 times and will use roots of unity to calculate the output.
    for i in range(N // 2):
        #   So as we are using global array for Roots which are generated previously
        # so, we use this equation to get correct nth root of unity in variable "RootIndex".
        RootIndex = i * (len(Roots) // (N))

        # Here we are getting our complex Nth Root real part in "Root" and 
        # Imaginary part in "Rooti".
        Root = Roots[RootIndex][0]
        # Here We are Taking Imaginary part as negative if We Want to perform IFFT.
        # because in IFFT we do Cos(Theta) - Sin(Theta).
        Rooti = Roots[RootIndex][1] if not IFFT else -Roots[RootIndex][1]

        # as per algo reference first we calculate Root * Odd output
        RootOdd, RootOddi = multiply_complex_num(Root, Rooti, Eodd[i][0], Eodd[i][1])

        # After which we Add or Subtract Even output base on index with Our Root * Odd.
        EvalReturn[i] = add_complex_num(Eeven[i][0], Eeven[i][1], RootOdd, RootOddi)
        EvalReturn[i + (N // 2)] = add_complex_num(
            Eeven[i][0], Eeven[i][1], -RootOdd, -RootOddi
        )

    # If We are doing IFFT We will divide our whole Output by 2.
    # we don't have divide function for Complex/tuples so we will multiply it by 1/2 = 0.5
    if IFFT:
        for i in range(N):
            EvalReturn[i] = multiply_complex_num(
                EvalReturn[i][0], EvalReturn[i][1], 0.5, 0
            )

    # we return merged Output array.
    return EvalReturn


def product_polynomial_evaluations(poly_1_evaluations, poly_2_evaluations):
    """This Function we use in order to multiply to point value representation in O(n) time.
    Here both array will be of same type due to previous setup/assumption of this code.

    Args:
        poly_1_evaluations (Array(Tuple[2])): point value representation for first polynomial.
        poly_2_evaluations (Array(Tuple[2])): point value representation for second polynomial.

    Returns:
        Array(Tuple[2]) : return multiplication of 2 point value representation.
    """
    multPoly = []
    for i in range(len(poly_1_evaluations)):
        multPoly.append(
            multiply_complex_num(
                poly_1_evaluations[i][0],
                poly_1_evaluations[i][1],
                poly_2_evaluations[i][0],
                poly_2_evaluations[i][1],
            )
        )
    return multPoly


def normalized(ComplexArray):
    """This Function we use for removing very small terms( like e^-18 )from our array and
    also to convert from complex number to real number by ignoring very small complex part.

    Args:
        ComplexArray (Array(tuple[2])): It's input array of complex number.

    Returns:
        array[int]: it return array of int.
    """
    for i in range(len(ComplexArray)):
        if ComplexArray[i][0] < 0.01 and ComplexArray[i][0] > -0.01:
            ComplexArray[i] = 0
        else:
            ComplexArray[i] = round(ComplexArray[i][0])
    return ComplexArray


if __name__ == "__main__":
    """
    This is our Main.
    """

    # Input for Degree of First Polynomial
    print("-" * 100)
    deg_1 = int(input("Please enter degree of Polynomial 1 : "))
    # Here we take space separated coefficient as input, remove any extra spaces at end 
    # then split it by space and map it to int after which we store it in array of integers.
    poly_1_in = list(
        map(
            int,
            input(
                "Please enter "
                + str(deg_1 + 1)
                + " coefficients for polynomial 1 in the increasing order of the degree of the monomials\n"
            )
            .strip()
            .split(),
        )
    )

    # Input for Degree of Second Polynomial
    print("-" * 100)
    deg_2 = int(input("Please enter degree of Polynomial 2 : "))
    # Here we take space separated coefficient as input, remove any extra spaces at end 
    # then split it by space and map it to int after which we store it in array of integers.
    poly_2_in = list(
        map(
            int,
            input(
                "Please enter "
                + str(deg_2 + 1)
                + " coefficients for polynomial 2 in the increasing order of the degree of the monomials\n"
            )
            .strip()
            .split(),
        )
    )

    # Here we are recalculating degree to avoid issues due to larger input.
    deg_1 = len(poly_1)
    deg_2 = len(poly_2)

    # Here we are initializing our input array with 0 and size is
    # grater than deg1 + deg2 + 1 which is in power of 2.
    # we do this because we have to evaluate our polynomial at minimum deg1 + deg2 + 1 points.
    poly_1 = [0] * find_N(deg_1, deg_2)
    poly_2 = [0] * find_N(deg_1, deg_2)

    # Here we are putting our Int Input to our defined Input array for each polynomial.
    for i in range(len(poly_1_in)):
        poly_1[i] = poly_1_in[i]

    for i in range(len(poly_2_in)):
        poly_2[i] = poly_2_in[i]

    # Here we are printing 2 polynomials.
    print("-" * 100)
    print("polynomial 1 : ", end="")
    print_poly(poly_1)
    print("polynomial 2 : ", end="")
    print_poly(poly_2)

    # Here we are calling for Naive Multiplication.
    naive_polynomial_multiplication(poly_1, poly_2)

    # Here we are printing output of Naive Multiplication stored in Global array naive_prod.
    print("-" * 100)
    print("polynomial Product (Naive) : ", end="")
    print_poly(naive_prod)

    # Here we are Generating N nth Roots of Unity In Global Array Roots.
    # We want to find number of Roots in power of 2 so we are using Function find_N for that.
    find_complex_roots(find_N(deg_1, deg_2))

    # Here we are finding Point Val representation for Both Polynomials.
    poly_1_evaluations = Eval(poly_1, IFFT=False)
    poly_2_evaluations = Eval(poly_2, IFFT=False)

    # Here we are performing Multiplication of 2 point value Representation.
    product_poly_evaluations = product_polynomial_evaluations(
        poly_1_evaluations, poly_2_evaluations
    )

    # Here we are Performing IFFT on Multiplication of two polynomials.
    # We are using normalized function on output to convert it to Array
    # of integers so we can print our resulting polynomial
    fft_prod = normalized(Eval(product_poly_evaluations, IFFT=True))

    # Here we are printing our resulting polynomial multiplied using FFT in O(nLog(n)) time.
    print("-" * 100)
    print("polynomial Product ( FFT ) : ", end="")
    print_poly(fft_prod)
    print("-" * 100)

    # This part of code is just to find out that both NAIVE and FFT are Outputting same result
    # for Testing purpose.
    for i in range(len(naive_prod)):
        if naive_prod[i] - fft_prod[i]:
            print("Both Naive And FFT Output Are Not Same, Maybe Some Error!")
            break
    else:
        print("Both Naive And FFT Output Are Same, ALL GOOD!")
    print("-" * 100)


""" Sample Output
----------------------------------------------------------------------------------------------------
Please enter degree of Polynomial 1 : 4
Please enter 5 coefficients for polynomial 1 in the increasing order of the degree of the monomials
5 4 8 7 9
----------------------------------------------------------------------------------------------------
Please enter degree of Polynomial 2 : 3
Please enter 4 coefficients for polynomial 2 in the increasing order of the degree of the monomials
2 5 4 8
----------------------------------------------------------------------------------------------------
polynomial 1 : 9x*4 + 7x*3 + 8x*2 + 4x + 5
polynomial 2 : 8x*3 + 4x*2 + 5x + 2
----------------------------------------------------------------------------------------------------
polynomial Product (Naive) : 72x*7 + 92x*6 + 137x*5 + 117x*4 + 110x*3 + 56x*2 + 33x + 10
----------------------------------------------------------------------------------------------------
polynomial Product ( FFT ) : 72x*7 + 92x*6 + 137x*5 + 117x*4 + 110x*3 + 56x*2 + 33x + 10
----------------------------------------------------------------------------------------------------
Both Naive And FFT Output Are Same, All Good!
----------------------------------------------------------------------------------------------------
"""
