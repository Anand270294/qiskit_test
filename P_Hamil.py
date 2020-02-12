from sympy import *

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
        self.obj_exp = sympify(self.expression_str)

        for i in range(len(var)):
            s_term = (1/2)*(I - self.Z[i])
            self.Hamil_exp = self.obj_exp.subs(self.variables[i],s_term)

        self.Hamil_exp = expand(self.Hamil_exp)
        self.Hamil_exp = self.Hamil_exp.subs(I,1)
        coeff = self.Hamil_exp.as_coefficients_dict()

        gbl_phse = coeff.get(1)
        self.Hamil_exp = self.Hamil_exp - gbl_phse

        # Convert to expression into a sympy poly to get list of monomial expression to build the QC
        self.quanCir_list = Poly(self.Hamil_exp).monoms()

    
    def get_ObjFun(self):
        return self.obj_exp

    def get_PHamil(self):
        return self.Hamil_exp

    def get_QClist(self):
        return self.quanCir_list
    


