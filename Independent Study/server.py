import accum
import ecc
import rsa
import socket
import sys
import threading

class Server:
    def __init__(self, accumulator, host='', port=50001, size=1024):
        self._accumulator = accumulator
        self._host = host
        self._port = port
        self._size = size
        self._acceptThread = threading.Thread(target=self.accept)
        self._values = { }

    def start(self):
        self._acceptThread.start()

    def accept(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.bind((self._host, self._port))
        s.listen(0)
        while True:
            client, address = s.accept()
            client.setsockopt(socket.IPPROTO_TCP, socket.TCP_NODELAY, 1)
            threading.Thread(target=self.echo, args=[client, address[0]]).start()

    def echo(self, client, address):
        print("Received connection from address:", address)
        data = client.recv(self._size).decode("utf-8").strip()
        print("Received message:", data)
        tokens = data.split()
        if len(tokens) > 0:
            command = tokens[0]
            if command == 'accumulate':
                value = eval(tokens[1])
                witness = self._accumulator.accumulate(value)
                for key in self._values:
                    self._values[address] = self._accumulator.update(self._values[address], value)
                self._values[address] = witness
                client.send((str(witness) + "\n").encode("utf-8"))
            elif command == 'verify':
                element, witness = eval(tokens[1]), eval(tokens[2])
                valid = self._accumulator.verify(witness, element)
                print(valid)
                client.send((str(valid) + "\n").encode("utf-8"))
        client.close()
        print("Closed connection to:", address)

if __name__ == '__main__':
    A = accum.Accumulator(rsa.g, rsa.exp, rsa.ver, rsa.update, [rsa.n])
    S = Server(A)
    S.start()
