// File:        crypto.h
// Author:      Liam Morris
// Description: Contains various supporting functions to be used in other
//              cryptography programs.

#include <map>
#include <vector>

// ##################
// DISCRETE LOGARITHM
// ##################
// Solves the discrete logarithm problem by brute force.
// Input:
//   int a - base for the exponentiation
//   int b - result of exponentiation
//   int p - modulo that we are working in
// Output:
//   int - the exponent that a is raised to that yields b
int discrete_log_bf(int a, int b, int p);

// Solves a portion of the discrete logarithm problem by the Pohlig-Hellman
// algorithm. Given a prime factor (q) and it's power (c), the result of the
// Chinese Remainder Theorem relevant to the prime factor will be returned.
// Input:
//   int a - base for the exponentiation
//   int b - result of exponentiation
//   int p - modulo that we are working in
//   int q - prime factor
//   int c - degree of prime factor
// Output:
//   int - result of Chinese Remainder Theorem relevant to prime factor
// TODO: Implement prime factorization function so that function can be called
//       with just a, b, and p as inputs (like other discrete log functions)
int pohlighellman(int a, int b, int p, int q, int c);

// Solves the discrete logarithm problem with Shank's algorithm.
// Input:
//   int a - base for the exponentiation
//   int b - result of exponentiation
//   int p - modulus that we are working in
// Output:
//   int - the exponent that a is raised to that yields b
int shanks(int a, int b, int p);

// ##################
// MODULAR ARITHMETIC
// ##################
// Finds the modular inverse of a number using the Extended Euclidean Algorithm.
// Input:
//   int a - the value for which the inverse is being found
//   int m - the modulus we are working in
// Output:
//   int - the modular inverse of a
int EEA(int a, int m);

// Solves the equation a^b mod p with the square-and-multiply algorithm.
// Input:
//   int a - base for exponentiation
//   int b - exponent
//   int p - modulus
// Output:
//   int - the result of the equation a^b mod p
int mod(int a, int b, int p);

// ######
// GROUPS
// ######
// Given a generator and a group, computes the order of the given generator.
// Input:
//   int g - the generator being used
//   int p - the modulo we are working in
// Output:
//   int - the order of the generator
int order(int g, int p);

// Given a group modulo p, computes the orders of all elements in the group.
// Input:
//   int p - the modulo we are working in
// Output:
//   map<int, int> - a map of <element, order> pairs in the group modulo p
std::map<int, int> compute_orders(int p);

// Generates a list of all powers of a given generator in a group modulo p.
// Input:
//   int a - the generator
//   int p - the modulo we are working in
// Output:
//   vector<int> - the different powers of the generator in the given group
std::vector<int> powers(int a, int p);

// ################
// FIELD ARITHMETIC
// ################
// Multiplies two polynomials together modulo x^2 + 1 in GF(p^2)
// Input:
//   std::vector<int> p1 - the first polynomial (ax + b)
//   std::vector<int> p2 - the second polynomial (cx + d)
//   int p - the prime used to generate the field
// Output:
//   std::vector<int> - the resulting polynomial
std::vector<int> gf_mul(std::vector<int> p1, std::vector<int> p2, int p);

// Raises a given polynomial to a given power modulo x^2 + 1 in GF(p^2)
// Input:
//   std::vector<int> p1 - the polynomial being raised to a power
//   int e - the power to which p1 is being raised
//   int p - the prime used to generate the field
// Output:
//   std::vector<int> - the resulting polynomial
std::vector<int> gf_exp(std::vector<int> p1, int e, int p);

// Computes the order of a given polynomial modulo x^2 + 1 in GF(p^2)
// Input:
//   std::vector<int> p1 - the polynomial whose order is being computed
//   int p - the prime used to generate the field
// Output:
//   int - the order of p1 modulo x^2 + 1 in GF(p^2)
int gf_order(std::vector<int> p1, int p);


// Computes the discrete logarithm of a given polynomial modulo x^2 + 1
// in GF(p^2)
// Input:
//   std::vector<int> a - the polynomial being used as a base
//   std::vector<int> b - the argument of the logarithm
//   int p - the generating prime for the finite field
// Output
//   int - the power to which a is raised such that b is returned
int gf_shanks(std::vector<int> a, std::vector<int> b, int p);

// Computes the multiplicative modular inverse of a given polynomial modulo
// x^2 + 1 in GF(p^2)
// Input:
//   std::vector<int> a - the polynomial whose inverse is being calculated
//   int p - the generating prime for the finite field
// Output:
//   std::vector<int> - the multiplicative inverse of a
std::vector<int> gf_inv(std::vector<int> a, int p);

// ###############
// ELLIPTIC CURVES
// ###############
// Represents an elliptic curve point.
struct point {
  // The coordinates of the elliptic curve point.
  int x, y;
  point() : x(0), y(0) { }
};

// Determines if a given value is a quadratic residue of a group.
// Input:
//   int y - the value being tested
//   int p - the modulus of the group
//   int (*exp)(int) - a function that represents the group's elliptic curve
// Output:
//   bool - whether or not the value is a quadratic residue
bool is_residue(int y, int p, int (*exp)(int));

// Computes the values for which there is a quadratic residue in a group.
// Input:
//   int p - the modulus of the group
//   int (*exp)(int) - a function that represents the group's elliptic curve
// Output:
//   vector<int> - the list of residues in the given group
std::vector<int> find_residues(int p, int (*exp)(int));

// Computes point addition of two elliptic curve points.
// Input:
//   point P - the first point
//   point Q - the second point
//   int p - the modulus of the group
//   int a - the value 'a' in the elliptic curve function x^2 + ax + b
// Output:
//   point - the resulting point of P + Q
point ec_add(point P, point Q, int p, int a = 1);

// Computes point multiplication of an elliptic curve point.
// Input:
//   point P - the point
//   int m - the number that the point is being multiplied by
//   int p - the modulus of the group
//   int a - the value 'a' in the elliptic curve function x^2 + ax + b
// Output:
//   point - the resulting point of mP
point ec_mul(point P, int m, int p, int a = 1);

// Converts an integer to it's non-adjacent form.
// Input:
//   int n - the number being coverted
// Output:
//   vector<int> - the NAF of the integer
std::vector<int> to_naf(int n);

// Computes point multiplication of an elliptic curve point using the DoubleAdd
// method.
// Input:
//   point P - the point
//   point m - the NAF of the number that the point is being multiplied by
//   int p - the modulus of the group
//   int a - the value 'a' in the elliptic curve function x^2 + ax + b
// Output:
//   point - the resulting point of mP
point ec_mul_dbl(point P, int m, int p, int a = 1);

// Computes the order of a point in its group.
// Input:
//   point P - the point
//   int p - the modulus of the group
//   int a - the value 'a' in the elliptic curve function x^2 + ax + b
// Output:
//   int - the order of point P in its group
int ec_order(point P, int p, int a = 1);

// Computes the square root of a value in its group.
// Input:
//   int a - the value for which the square root is being computed.
//   int p - the modulus of the group where p = 3 mod 4
// Output:
//   int - the square root of the value
int ec_sqrt(int a, int p);

// Computes all possible coordinates of an elliptic curve group.
// Input:
//   int (*exp)(int) - the function that represents the elliptic curve
//   int p - the modulus of the group
// Output:
//   vector<point> - a list of points on the curve
std::vector<point> coordinates(int (*exp)(int), int p);

// Compresses an elliptic curve point.
// Input:
//   point P - the point being compressed
// Output:
//   point - the compressed form of P
point ec_compress(point P);

// Decompresses an elliptic curve point.
// Input:
//   point P - the point being decompressed
//   int p - the modulus of the group
//   int (*exp)(int) - the function that represents the elliptic curve
// Output:
//   point the decompressed form of P
point ec_decompress(point P, int p, int (*exp)(int));
