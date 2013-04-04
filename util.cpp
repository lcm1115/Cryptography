// File:        util.cpp
// Author:      Liam Morris
// Description: Implements the functions prototyped in util.h.

#include "util.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stdio.h>
#include <utility>
#include <vector>

using namespace std;

// ##################
// DISCRETE LOGARITHM
// ##################
int discrete_log_bf(int a, int b, int p) {
  int val = 1;
  // For all values from 0 to p-1, see if the discrete log is found.
  for (int i = 0; i < p; ++i) {
    if (val == b) return i;
    val *= a;
    val %= p;
  }
}

int pohlighellman(int a, int b, int p, int* q_vals, int* c_vals, int num_vals, bool is_squared = false) {
  int result = 0;
  int inv = EEA(a, p);
  int c, exp, i, m, mul, q, temp_result;
  int pval = p;
  if (is_squared) {
    pval *= p;
    pval -= p;
    p = pow(p, 2);
  }
  for (int k = 0; k < num_vals; ++k) {
    int B = b;
    temp_result = 0;
    // Get current q and c values
    q = q_vals[k];
    c = c_vals[k];
    printf("Computing with q = %d and c = %d\n", q, c);
    for (int j = 0; j < c; ++j) {
      printf("j = %d\tB = %d\n", j, B);
      exp = pval / pow(q, j + 1);
      m = mod(B, exp, p);
      printf("\tm = %d", m);
      // Find i value
      for (i = 0; i < p; ++i) {
        if (mod(a, i * pval / q, p) == m) {
          break;
        }
      }
      i %= q;
      printf("\ti = %d\n", i);
      // Add value to temp_result and get new B value
      temp_result += i * pow(q, j);
      mul = mod(inv, i * pow(q, j), p);
      B = (B * mul) % p;
    }
    // Add result of this q value to end result
    m = (pval) / pow(q, c);
    temp_result *= m * EEA(m, pow(q, c));
    temp_result %= pval;
    printf("Result from CRT to be added to total: %d\n\n", temp_result);
    result += temp_result;
    result %= pval;
  }

  return result;
}

int shanks(int a, int b, int p) {
  // Compute m value
  int m = ceil(sqrt(p - 1));

  // Compute a^m % p
  int mul = mod(a, m, p);

  // Populate L1 (as a map)
  int val = 1;
  map<int, int> L1;
  for (int i = 0; i < m; ++i) {
    L1[val] = i;
    val *= mul;
    val %= p;
  }

  // Reset val
  val = b;

  // Compute values of L2 until one is found in L!
  int a1 = 1;
  for (int i = 0; i < m; ++i) {
    // If current value is found in L1, return the answer (mj + i)
    if (L1.find(val) != L1.end()) {
      return L1[val] * m + i;
    }
    a1 *= a;
    a1 %= p;
    val = b * EEA(a1, p);
    val %= p;
  }

  return -1;
}

// ##################
// MODULAR ARITHMETIC
// ##################
int EEA(int a, int m) {
  int x = m, y = a, t2 = 0, t1 = 1;
  while (true) {
    int r = x % y;
    if (r == 0) break;
    int q = (x - r) / y;
    int t = t2 - t1 * q;
    x = y, y = r, t2 = t1, t1 = t;
  }
  if (y != 1) return -1;
  if (t1 < 0) {
    t1 += m;
  }
  return t1;
}

int mod(int a, int b, int p) {
  int mod = 1;
  int mul = a;
  while (b > 0) {
    if (b & 1) {
      mod *= mul;
      mod %= p;
    }
    mul *= mul;
    mul %= p;
    b = b >> 1;
  }
  return mod;
}

int order(int g, int p) {
  int val = (g * g) % p;
  int order = 1;
  while (val != g && val > 0) {
    order++;
    val *= g;
    val %= p;
  }

  return order;
}

map<int, int> compute_orders(int p) {
  map<int, int> orders;
  for (int i = 2; i < p; ++i) {
    orders[i] = order(i, p);
  }
  return orders;
}

vector<int> powers(int a, int p) {
  vector<int> powers;
  int val = 1;
  powers.push_back(val);
  val = a;
  while (val != 1) {
    powers.push_back(val);
    val *= a;
    val %= p;
  }

  return powers;
}

// ################
// FIELD ARITHMETIC
// ################
vector<int> gf_mul(vector<int> p1, vector<int> p2, int p) {
  vector<int> result(2);
  result.at(0) = abs(p1.at(0) * p2.at(1) + p1.at(1) * p2.at(0)) % p;
  result.at(1) = abs(p1.at(1) * p2.at(1) - p1.at(0) * p2.at(0)) % p;
  return result;
}

vector<int> gf_exp(vector<int> p1, int e, int p) {
  vector<int> result;
  result.push_back(0);
  result.push_back(1);
  for (int i = 0; i < e; ++i) {
    result = gf_mul(result, p1, p);
  }
  return result;
}

int gf_order(vector<int> p1, int p) {
  int order = 1;
  vector<int> cur(p1);
  while ((cur.at(0) != 0 || cur.at(1) != 1) && order < p * p) {
    cur = gf_mul(cur, p1, p);
    ++order;
  }
  return order;
}
