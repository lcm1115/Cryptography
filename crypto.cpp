// File:        crypto.cpp
// Author:      Liam Morris
// Description: Implements the functions prototyped in crypto.h.

#include "crypto.h"

#include <algorithm>
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

    return -1;
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

    // Compute values of L2 until one is found in L1
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
    if (a == 0) return 0;
    if (a < 0) a += m;
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

    // Iterate across all bits in the number.
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
    int a = (p1.at(0) * p2.at(1) + p1.at(1) * p2.at(0)) % p;
    int b = (p1.at(1) * p2.at(1) - p1.at(0) * p2.at(0)) % p;
    while (a < 0) {
        a += p;
    }
    while (b < 0) {
        b += p;
    }
    result.at(0) = a;
    result.at(1) = b;
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
        if (cur.at(0) == p1.at(0) && cur.at(1) == p1.at(1)) {
            break;
        }
        ++order;
    }
    return order;
}

int gf_shanks(vector<int> a, vector<int> b, int p) {
    int m = p;
    vector<int> val(2);
    val.at(0) = 0;
    val.at(1) = 1;
    vector<vector<int> > L1;
    for (int i = 0; i < m; ++i) {
        L1.push_back(val);
        val = gf_mul(val, a, p); 
    }
    vector<int> inv = gf_exp(gf_inv(a, p), m, p);
    vector<int> y = b;
    for(int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            if (L1[j] == y) {
                return (i * m + j);
            }
        }
        y = gf_mul(y, inv, p);
    }

    return -1;
}

vector<int> gf_inv(vector<int> a, int p) {
    vector<int> poly(2);
    vector<int> result;
    for (int i = 0; i < p; ++i) {
        for (int j = 0; j < p; ++j) {
            poly.at(0) = i;
            poly.at(1) = j;
            result = gf_mul(a, poly, p);
            if (result.at(0) == 0 && result.at(1) == 1) {
                return poly;
            }
        }
    }

    return a;
}

// ###############
// ELLIPTIC CURVES
// ###############
bool is_residue(int y, int p, int (*exp)(int)) {
    // Euler criterion for quadratic residue
    return mod(y, (p - 1) / 2, p) == 1;
}

vector<int> find_residues(int p, int (*exp)(int)) {
    // Iterate across all values in the group and determine which are residues
    vector<int> residues;
    for (int i = 1; i < p; ++i) {
        int y = (*exp)(i);
        // If a value is a residue, add it to the list
        if (is_residue(y, p, exp)) {
            residues.push_back(i);
        }
    }

    return residues;
}

point ec_add(point P, point Q, int p, int a) {
    // If either point is 0, return the other point
    if (P.x == 0 && P.y == 0) return Q;
    if (Q.x == 0 && Q.y == 0) return P;

    point result;
    result.x = 0;
    result.y = 0;

    // If the points are inverses of each other, return 0
    if (P.x == Q.x &&
            (P.y + Q.y == p || (P.y == 0 && Q.y == 0))) return result;
    int lambda;
    int t;

    // Compute lambda depending on P and Q equality.
    if (P.x != Q.x || P.y != Q.y) {
        lambda = ((Q.y - P.y) * EEA(Q.x - P.x, p)) % p;
    } else {
        lambda = ((3 * P.x * P.x + a) * EEA(2 * P.y, p)) % p;
    }

    // Compute the x and y values of the resulting point
    result.x = (lambda * lambda - P.x - Q.x) % p;
    result.y = (lambda * (P.x - result.x) - P.y) % p;

    // Correct negative values
    if (result.x < 0) result.x += p;
    if (result.y < 0) result.y += p;

    return result;
}

point ec_mul(point P, int m, int p, int a) {
    point result;
    result.x = 0;
    result.y = 0;

    // Add point m times
    for (int i = 0; i < m; ++i) {
        result = ec_add(P, result, p, a);
    }

    return result;
}

vector<int> to_naf(int n) {
    // Convert integer to standard binary form
    vector<int> digits;
    while (n > 0) {
        digits.push_back(n & 1);
        n = n >> 1;
    }

    // Iterate across digits and convert as needed
    for (int i = 0; i < digits.size() - 1; ++i) {
        int j = i + 1;

        // Identify block of 1 digits
        if (digits.at(i) == 1 && digits.at(j) == 1) {
            digits.at(i) = -1;

            // Convert all 1's to 0's until block is in NAF
            while (digits.at(j) == 1) {
                digits.at(j) = 0;
                ++j;

                // If another digit is needed, add it
                if (j == digits.size()) {
                    digits.push_back(0);
                }
            }

            // Convert last digit in the block to 1
            digits.at(j) = 1;
        }
    }
    return digits;
}

point ec_mul_dbl(point P, int m, int p, int a) {
    point Q;
    point P_neg(P);
    P_neg.y = p - P_neg.y;

    // Get NAF of multiplier
    vector<int> NAF = to_naf(m);

    for (int i = NAF.size() - 1; i >= 0; --i) {
        // Double the point
        Q = ec_add(Q, Q, p);

        // If the current index is a 1, add the point
        if (NAF.at(i) == 1) {
            Q = ec_add(Q, P, p);
        }
        // If the current index is a -1, subtract the point
        else if (NAF.at(i) == -1) {
            Q = ec_add(Q, P_neg, p);
        }
    }

    return Q;
}

int ec_order(point P, int p, int a) {
    // Initialize order to 1 and point to 0
    int order = 1;
    point result;
    result.x = P.x;
    result.y = P.y;

    // Add the point to the current point until 0 is reached
    while (result.x != 0 || result.y != 0) {
        result = ec_add(P, result, p, a);
        ++order;
    }

    return order;
}

int ec_sqrt(int a, int p) {
    // Use formula to compute square root
    return mod(a, (p + 1) / 4, p);
}

vector<point> coordinates(int (*exp)(int), int p) {
    vector<point> coordinates;
    point curr;

    // Add 0 to coordinates
    coordinates.push_back(curr);
    for (int i = 0; i < p; ++i) {
        int y1 = (*exp)(i);
        int m = mod(y1, (p - 1) / 2, p);

        // If the value is a quadratic residue, compute both points and add them to
        // the list
        if (m == 1) {
            y1 = ec_sqrt(y1, p);
            int y2 = p - y1;
            if (y2 < 0) y2 += p;
            curr.x = i;
            curr.y = min(y1, y2);
            coordinates.push_back(curr);
            curr.y = max(y1, y2);
            coordinates.push_back(curr);
        }
        // If the resulting value is 0, add the one corresponding point to the list
        else if (m == 0) {
            curr.x = i;
            curr.y = 0;
            coordinates.push_back(curr);
        }
    }

    return coordinates;
}

point ec_compress(point P) {
    point result(P);
    result.y %= 2;
    return result;
}

point ec_decompress(point P, int p, int (*exp)(int)) {
    int x = P.x;
    int z = (*exp)(x);
    if (is_residue(z, p, exp)) {
        int y = ec_sqrt(z, p);
        if (y % 2 == P.y) {
            P.y = y;
            return P;
        }
        else {
            P.y = p - y;
            return P;
        }
    }

    return P;
}
