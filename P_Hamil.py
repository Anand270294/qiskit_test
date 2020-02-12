from sympy import *
from sympy.core.compatibility import exec_
from sympy.parsing.sympy_parser import parse_expr

class P_Hamil:
    """ Class for converting Objective functions into Phase/Cost Hamiltonians,"""
    
    

    def __init__(self, expr_str, var):
        self.__checkVariables(expr_str, var)
        
        # Initialize expressions and variables
        self.expression_str = expr_str
        self.variables = var



    def __checkVariables(self, expression, variables):
        for v in variables:
            if(v in expression):
                pass
            else:
                raise ValueError('Variables Mismatch! Unable to find {} in the Objective Function: {}'.format(v, expression))

    
    def Hamify(self):
        self.H_exp = parse_expr(self.expression_str)

        print(self.H_exp)


