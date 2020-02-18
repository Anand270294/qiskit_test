from sympy import *
from qiskit import Aer, IBMQ
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit,execute

class M_Hamil:
    """ Class for generating Mixer Hamiltonians for the phase/cost Hamiltonians """
    I = symbols('I')
    X = []
    Y = []
    Z = [] 

    obj_exp = ""
    variables = []

    def __init__(self, expr_str, var):
        self.obj_exp = sympify(expr_str)
        self.variables = var

        for i in range(len(var)):
            self.Z.append(symbols('Z_{}'.format(var[i][len(var[i])-1])))
            self.Y.append(symbols('Y_{}'.format(var[i][len(var[i])-1])))
            self.X.append(symbols('Y_{}'.format(var[i][len(var[i])-1])))

    
    # Mixer Hamiltonian for qubits to have dynamicity between {0,1}
    def generalXMixer(self, beta, q, measure=False):
        self.mixer_circuit = QuantumCircuit(len(self.variables), len(self.variables))
        
        for i in range(len(self.variables)):
            self.mixer_circuit.rx(beta, i)

        if measure == True:
            self.mixer_circuit.barrier()
            self.mixer_circuit.measure(range(len(self.variables)), range(len(self.variables)))


        return self.mixer_circuit


