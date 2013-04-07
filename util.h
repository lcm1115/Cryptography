// File:        util.h
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
//   vector<int> p1 - the first polynomial (ax + b)
//   vector<int> p2 - the second polynomial (cx + d)
//   int p - the prime used to generate the field
// Output:
//   vector<int> - the resulting polynomial
std::vector<int> gf_mul(std::vector<int> p1, std::vector<int> p2, int p);

// Raises a given polynomial to a given power modulo x^2 + 1 in GF(p^2)
// Input:
//   vector<int> p1 - the polynomial being raised to a power
//   int e - the power to which p1 is being raised
//   int p - the prime used to generate the field
// Output:
//   vector<int> - the resulting polynomial
std::vector<int> gf_exp(std::vector<int> p1, int e, int p);

// Computes the order of a given polynomial modulo x^2 + 1 in GF(p^2)
// Input:
//   vector<int> p1 - the polynomial whose order is being computed
//   int p - the prime used to generate the field
// Output:
//   int - the order of p1 modulo x^2 + 1 in GF(p^2)
int gf_order(std::vector<int> p1, int p);

int gf_shanks(std::vector<int> a, std::vector<int> b, int p);

std::vector<int> gf_inv(std::vector<int> a, int p);
