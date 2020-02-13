from sympy import *
from qiskit import Aer, IBMQ
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit,execute

class M_Hamil:
    """ Class for generating Mixer Hamiltonians for the phase/cost Hamiltonians """
    