def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)

def modinv(a, m):
    a %= m
    g, x, y = egcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m

def add(p, q, a, n):
    if p == (0, 0):
        return q
    elif q == (0, 0):
        return p
    else:
        if p[0] != q[0]:
            L = (p[1] - q[1]) * (modinv(p[0] - q[0], n))
        else:
            if p[1] != q[1] or q[1] == 0:
                return (0, 0)
            L = ((3 * q[0] ** 2 + a) * modinv(2 * q[1], n)) % n
        x = (L ** 2 - p[0] - q[0]) % n
        y = ((q[0] - x) * L - q[1]) % n
        return (x, y)

def mul(p, k, a, n):
    q = (0, 0)
    for i in range(k):
        q = add(p, q, a, n)
    return q

def acc(p, k, data):
    return mul(p, k, data[0], data[1])

def ver(value, witness, element, data):
    return value == mul(witness, element, data[0], data[1])

def update(witness, element, data):
    return mul(witness, element, data[0], data[1])
