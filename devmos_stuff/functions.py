"""Provide the primary functions."""

import numpy as np
import math    
import networkx as nx
import matplotlib.pyplot as plt
import random
import copy as cp         
random.seed(2)

class BitString:
    def __init__(self, N):
        self.N = N
        self.config = np.zeros(N, dtype=int) 

    def __repr__(self):
        out = ""
        for i in self.config:
            out += str(i)
        return out

    def __eq__(self, other):        
        return all(self.config == other.config)
    
    def __len__(self):
        return len(self.config)

    def on(self):
        value = 0
        for char in self.config:
            if char == "1" or char == 1:
                value += 1
        return value
        """
        Return number of bits that are on
        """

    def off(self):
        value = 0
        for char in self.config:
            if char == "0" or char == 0:
                value += 1
        return value
        """
        Return number of bits that are on
        """

    def flip_site(self,i):
        if self.config[i] == 1:
            self.config[i] = 0
        else:
            self.config[i] = 1
        """
        Flip the bit at site i
        """
    
    def integer(self):
        iter = 0
        exp = 1
        while exp <= len(self.config):
            if self.config[-exp] != None:
                if self.config[-exp] == 1 or self.config[-exp] == "1":
                    iter += 2 ** (exp-1)
                    exp += 1
                else:
                    exp += 1
            else:
                exp += 1
        return iter

        """
        Return the decimal integer corresponding to BitString
        """
 

    def set_config(self, s:list[int]):
        self.config = ""
        for value in s:
            self.config += str(value)
        """
        Set the config from a list of integers
        """

    def set_integer_config(self, dec:int, digits=None):
        n = 1
        if dec == 0 and digits == None:
            n = 1
        elif digits != None:
            n = digits
        else:
            n = math.ceil(math.log2(dec + 1))
        array = np.zeros(6, dtype=int)
        if dec == 0 or dec == None:
            softcap = 0
            dec = 0
        elif dec == 1:
            softcap = 1
        else:
            softcap = math.ceil(math.log(dec,2))
        if dec % 2 == 1 and dec != 0:
            array[-1] = 1
        while softcap > 0:
            if (dec / (2 ** softcap)) >= 1:
                array[-softcap-1] = 1
                dec -= 2 ** softcap
                softcap -= 1
            else:
                try:
                    array[-softcap-1] = 0
                    softcap -= 1
                except IndexError:
                    softcap -= 1
        try: 
            self.config = array
            return array
        except AttributeError:
            return array
        

class IsingHamiltonian:
    def __init__(self, G):
        self.G = nx.adjacency_matrix(G).todense()
        self.mu = []
    
    def energy(self, bs):
        bs = str(bs)
        Energy = 0
        row_index = 0
        for row in self.G:
            column_index = 0
            for value in row:
                if value != 0:
                    Spin1 = 0
                    Spin2 = 0
                    if bs[row_index] == "1":
                        Spin1 = value
                    elif bs[row_index] == "0":
                        Spin1 = -value
                    if bs[column_index] == "1":
                        Spin2 = value
                    elif bs[column_index] == "0":
                        Spin2 = -value 
                    Energy += ((Spin1 * Spin2) / (abs(Spin1) + abs(Spin2)))
                    column_index += 1
                else:
                    column_index += 1
            row_index += 1
        iter = 0
        for char in bs:
            if char == "1" and len(self.mu) != 0:
                Energy += self.mu[iter]
            elif char == "0" and len(self.mu) != 0:
                Energy -= self.mu[iter]
            else:
                pass
            iter += 1
        return Energy

    
    def compute_average_values(self, T:float=1, bs=BitString):

        E  = 0.0
        M  = 0.0
        E2EXP = 0.0
        M2EXP = 0.0
        Z = 0.0
        HC = 0.0
        MS = 0.0
        k = 1
        B = 1 / (k * T)
    
        # Calculating Z
        for index in range(64):
            if index == 0:
                indexcatch = 1
            else:
                indexcatch = index
            BitString.set_integer_config(bs, index)
            Z += math.e ** (-B * float(IsingHamiltonian.energy(self, bs)))

        # Calculating Primary Values
        for index in range(64):
            if index == 0:
                indexcatch = 1
            else:
                indexcatch = index
            BitString.set_integer_config(bs, index)
            E += (float(IsingHamiltonian.energy(self, bs)) * (math.e ** (-B * float(IsingHamiltonian.energy(self, bs)))))/Z
            E2EXP += (((float(IsingHamiltonian.energy(self, bs))) ** 2) * (math.e ** (-B * float(IsingHamiltonian.energy(self, bs)))))/Z
            M += ((BitString.on(self) - BitString.off(self)) * (math.e ** (-B * float(IsingHamiltonian.energy(self, bs)))))/Z
            M2EXP += (((BitString.on(self) - BitString.off(self)) ** 2) * (math.e ** (-B * float(IsingHamiltonian.energy(self, bs)))))/Z

        HC = (E2EXP - (E ** 2)) * (T ** -2)
        MS = (M2EXP - (M ** 2)) * (T ** -1)
    
        return E, M, HC, MS
    
    def set_mu(self, array):
        self.mu = array
    


class MonteCarlo:
    def __init__(self, ham, G):
        self.ham = ham
        self.G = nx.adjacency_matrix(G).todense()
        self.mu = []

    def set_mu(self, array):
        self.mu = array

    def run(self, T, turns, burns=0):
        E = []
        M = []
        LowestE = None
        tempE = None
        tempM = None
        while turns != 0:
            bs = BitString.set_integer_config(None, random.randint(0,63))
            tempE,tempM = IsingHamiltonian.compute_average_values(self, T, bs)
            if tempE <= LowestE or LowestE == None:
                    LowestE = tempE 
            if burns != 0:
                burns -= 1
                turns -= 1
                tempE = None
                pass
            else:
                E.append(tempE)
                M.append(tempM)
                tempE = None
                turns -= 1
            
    

def canvas(with_attribution=True):
    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())
