import math
import sys
import matplotlib.pyplot as mp
from itertools import *


class NS:
    def __init__(self):
        self._xran = []
        self._n = 0.3
        self._c = []
        for index in range(16):
            vector = bin(index)[2:].zfill(4)
            x1 = int(vector[0])
            x2 = int(vector[1])
            x3 = int(vector[2])
            x4 = int(vector[3])
            self._xran.append(self.func(x1, x2, x3, x4))
            if self._xran[-1] == 1:
                self._c.append([x1, x2, x3, x4])
        self._j = len(self._c)
        self._w = [0.0 for i in range(self._j + 1)]
        self._e = []

    @staticmethod
    def func(x1, x2, x3, x4):
        x1 = True if x1 == 1 else False
        x2 = True if x2 == 1 else False
        x3 = True if x3 == 1 else False
        x4 = True if x4 == 1 else False
        f = (x1 or x2) and x3 or x4
        return 1 if f is True else 0

    def Tru_ta(self):
        print(' x1 | x2 | x3 | x4 | F ')
        for index in range(16):
            vector = bin(index)[2:].zfill(4)
            x1 = int(vector[0])
            x2 = int(vector[1])
            x3 = int(vector[2])
            x4 = int(vector[3])
            print('  {} |  {} |  {} |  {} | {} '.format(x1, x2, x3, x4, self._xran[index]))

    def tra_net(self, net):
        return 1 if net > (- sys.float_info.epsilon) else 0

    def log_net(self, net):
        return 1/2*(math.tanh(net)+1)
    
    def Fi_g(self, x, n):
        f = 0
        for i in range(4):
            f += math.pow(x[i] - self._c[n][i], 2)
        return math.exp(-f)

    def comb(self, w, F):
        y = [0 for i in range(16)]
        for index in range(16):
            fi = []
            vector = bin(index)[2:].zfill(4)
            x1 = int(vector[0])
            x2 = int(vector[1])
            x3 = int(vector[2])
            x4 = int(vector[3])
            fi.append(1)
            for i in range(1, self._j + 1):
                fi.append(self.Fi_g([x1, x2, x3, x4], i - 1))
            net = w[0]
            for i in range(1, self._j + 1):
                net += w[i] * fi[i]
            if F == 1:
                y[index] = self.tra_net(net)
            elif F == 3:
                if self.log_net(net) >= 0.5:
                    y[index] = 1
                else:
                    y[index] = 0
            delta = self._xran[index] - y[index]
            if delta != 0:
                return False
        return True

    def edu(self, F):
        all_vectors = [i for i in range(16)]
        self.pr_edu(all_vectors, F)

    def pr_edu(self, vectors, F):
        self._e = []
        self._w = [0.0 for i in range(self._j + 1)]
        era = 0
        y = [0 for i in range(len(vectors))]
        line = ''
        if len(vectors) < 16:
            line += '{} vectors'.format(len(vectors))
        if F == 1:
            print('Threshold {}'.format(line).center(106))
        else:
            print('Logistic FA {}'.format(line).center(106))
        print('Era'.center(11), end='|')
        print('Weight'.center(62), end='                                              |')
        print('Out'.center(47), end='|')
        print('ERR'.center(11))
        while era < 200:
            print('{:2d}'.format(era).center(11), end='|')
            str1 = ''
            for i in range(len(self._w)):
                if F == 1:
                    str1 += '{:8.3} '.format(self._w[i])
                else:
                    str1 += '{:8.4} '.format(self._w[i])
            str1.center(56)
            print(str1, end='|')
            E = 0
            for index in range(len(vectors)):
                fi = []
                vector = bin(vectors[index])[2:].zfill(4)
                x1 = int(vector[0])
                x2 = int(vector[1])
                x3 = int(vector[2])
                x4 = int(vector[3])
                fi.append(1)
                for i in range(1, self._j + 1):
                    fi.append(self.Fi_g([x1, x2, x3, x4], i - 1))
                net = self._w[0]
                for i in range(1, self._j + 1):
                    net += self._w[i] * fi[i]
                if F == 1:
                    y[index] = self.tra_net(net)
                elif F == 3:
                    if self.log_net(net) >= 0.5:
                        y[index] = 1
                    else:
                        y[index] = 0
                d = self._xran[vectors[index]] - y[index]
                if d != 0:
                    E += 1
                    for i in range(self._j + 1):
                        dw = self._n * d * fi[i]
                        if F == 3:
                            dw *= self.log_net(net) * (1 - self.log_net(net))
                        self._w[i] += dw
            self._e.append(E)
            era += 1
            Vec_o = ''
            for j in range(len(vectors) - 1):
                Vec_o += ' {},'.format(y[j])
            Vec_o += ' {}'.format(y[-1])
            print(Vec_o.center(47), end='|')
            print('{}'.format(E).center(11))
            if E == 0:
                print()
                break

    def plot(self):
        eras = [i for i in range(len(self._e))]
        mp.xlabel('Era')
        mp.ylabel('ERR')
        mp.grid()
        mp.plot(eras, self._e)
        mp.show()

    def min_s(self, FA):
        comb = []
        a_vec = []
        for index in range(16):
            a_vec.append(index)
        for i in range(2, 16):
            for vecs in combinations(a_vec, i):
                w = [0.0 for index in range(self._j + 1)]
                y = [0 for index in range(i)]
                era = 0
                while era < 50:
                    E = 0
                    for index in range(len(vecs)):
                        fi = []
                        vector = bin(vecs[index])[2:].zfill(4)
                        x1 = int(vector[0])
                        x2 = int(vector[1])
                        x3 = int(vector[2])
                        x4 = int(vector[3])
                        fi.append(1)
                        for j in range(1, self._j + 1):
                            fi.append(self.Fi_g([x1, x2, x3, x4], j - 1))
                        net = w[0]
                        for j in range(1, self._j + 1):
                            net += w[j] * fi[j]
                        if FA == 1:
                            y[index] = self.tra_net(net)
                        elif FA == 3:
                            if self.log_net(net) >= 0.5:
                                y[index] = 1
                            else:
                                y[index] = 0
                        d = self._xran[vecs[index]] - y[index]
                        if d != 0:
                            E += 1
                            for j in range(self._j + 1):
                                dw = self._n * d * fi[j]
                                if FA == 3:
                                    dw *= self.log_net(net) * (1 - self.log_net(net))
                                w[j] += dw
                    era += 1
                    if E == 0:
                        podh = self.comb(w, FA)
                        if podh:
                            comb = vecs
                            break
                if len(comb) > 0:
                    break
            if len(comb) > 0:
                break
        print(comb, end=' = (')
        for vector in comb:
            bin_vector = bin(vector)[2:].zfill(4)
            print(' {} '.format(bin_vector), end='')
        print(')')
        self.pr_edu(comb, FA)



nei = NS()
nei.Tru_ta()
nei.edu(1)
nei.plot()
nei.min_s(1)
nei.plot()
nei2 = NS()
nei2.edu(3)
nei2.plot()
nei2.min_s(3)
nei2.plot()


