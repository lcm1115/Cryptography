import ecc
import rsa

class Accumulator:
    def __init__(self, value, accumulate, verify, update, data=[]):
        self._value = value
        self._accumulate = accumulate
        self._verify = verify
        self._update = update
        self._data = data

    def accumulate(self, element):
        witness = self._value
        self._value = self._accumulate(self._value, element, self._data)
        return witness

    def verify(self, witness, element):
        return self._verify(self._value, witness, element, self._data)

    def update(self, witness, element):
        return self._update(witness, element, self._data)
