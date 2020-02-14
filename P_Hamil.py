from sympy import *
from qiskit import Aer, IBMQ
from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit,execute

class P_Hamil:
    """ Class for converting Objective functions into Phase/Cost Hamiltonians,"""
    I = symbols('I')
    X = []
    Y = []
    Z = []    

    obj_exp = ""
    Hamil_exp = ""
    quanCir_list = []
    

    def __init__(self, expr_str, var):
        self.__checkVariables(expr_str, var)
        
        # Initialize expressions and variables
        self.expression_str = expr_str
        self.variables = var

        self.obj_exp = sympify(self.expression_str)

        # Number of pauli objects are limited to the number of variables/qubits being used
        for i in range(len(var)):
            self.Z.append(symbols('Z_{}'.format(var[i][len(var[i])-1])))
            self.Y.append(symbols('Y_{}'.format(var[i][len(var[i])-1])))
            self.X.append(symbols('Y_{}'.format(var[i][len(var[i])-1])))



    def __checkVariables(self, expression, variables):
        for v in variables:
            if(v in expression):
                pass
            else:
                raise ValueError('Variables Mismatch! Unable to find {} in the Objective Function: {}'.format(v, expression))

    
    def Hamify(self):
        self.Hamil_exp = self.obj_exp

        for i in range(len(self.variables)):
            s_term = (1/2)*(I - self.Z[i])
            self.Hamil_exp = self.Hamil_exp.subs(self.variables[i],s_term)

        self.Hamil_exp = expand(self.Hamil_exp)
        self.Hamil_exp = self.Hamil_exp.subs(I,1)
        coeff = self.Hamil_exp.as_coefficients_dict()

        gbl_phse = coeff.get(1)
        self.Hamil_exp = self.Hamil_exp - gbl_phse

        # Convert to expression into a sympy poly to get list of monomial expression to build the QC
        self.quanCir_list = Poly(self.Hamil_exp).monoms()

    

    # Map the qubits directly to each variable 
    def perQubitMap(self, gamma, p, barrier=False, initial_Hadamard=False):
        self.p_hamilCir = QuantumCircuit(len(self.variables),len(self.variables))

        if initial_Hadamard == True:
            for i in range(len(self.variables)):
                self.p_hamilCir.h(i)
            self.p_hamilCir.barrier()

        for sub_cir in self.quanCir_list:
            if sum(sub_cir) > 1:
                indices = [i for i, x in enumerate(sub_cir) if x == 1]
                for i in range(len(indices)):
                    if i == len(indices) - 1:
                        self.p_hamilCir.rz(gamma, indices[i])
                    else:
                        self.p_hamilCir.cx(indices[i], indices[i + 1])
                
                for i in range(len(indices)-1,0,-1):
                    self.p_hamilCir.cx(indices[i-1], indices[i])

            else:
                self.p_hamilCir.rz(gamma, sub_cir.index(1))


        if barrier == True:
            self.p_hamilCir.barrier()

        return self.p_hamilCir

    # def perEdgeMap(self):
    #     continue

    def expectation_value(self, results, shots):
        counts = results.get_counts()

        i = 0
        expectation = 0
        for c in counts:
            x = self.obj_exp
            for i in range(len(self.variables)):
                x = x.subs(self.variables[i], int(c[i]))
                
            expectation += (counts[c] / shots) * x
            
        return expectation

    def drawCircuit(self,output):
        return self.p_hamilCir.draw(output=output)

    def get_ObjFun(self):
        return self.obj_exp

    def get_PHamil(self):
        return self.Hamil_exp

    def get_QClist(self):
        return self.quanCir_list
    


