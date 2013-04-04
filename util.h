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
int order(int g, int p);

std::map<int, int> compute_orders(int p);

// ################
// FIELD ARITHMETIC
// ################
std::vector<int> powers(int a, int p);
std::vector<int> gf_mul(std::vector<int> p1, std::vector<int> p2, int p);
std::vector<int> gf_exp(std::vector<int> p1, int e, int p);
int gf_order(std::vector<int> p1, int p);
