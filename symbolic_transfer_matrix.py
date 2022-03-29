from sympy import re, im, I, symbols, Matrix, exp, eye, simplify
from numpy import pi
import sys

class TransferMatrix():

    def symPropagationMatrix(self, layer_index: int):

        layer_index_str = str(layer_index)
        phi = symbols('phi{sub}'.format(sub = layer_index_str))

        P = Matrix([[exp(I*phi), 0], 
                    [0, exp(-I*phi)]])

        return P

    def symInterfaceMatrix(self, interface_index: int):

        right_index_str = str(interface_index)
        left_index_str = str(interface_index - 1)

        r, t = symbols('r{sub0}{sub1} t{sub0}{sub1}'.format(sub0 = left_index_str, 
                                                            sub1 = right_index_str))

        D = Matrix([[1, r], 
                    [r, 1]])/t
        
        return D
    
    def symFullMatrix(self, num_layers: int):

        """
        Symbolic calculation of the Transfer Matrix for n layers.

        """
        M = eye(2)
        num_interfaces = num_layers - 1

        M_list = []  
        D_list = [] 
        P_list = [] 
        
        for i in reversed(range(1, num_interfaces + 1)):

            if i == 1:
               
                D = self.symInterfaceMatrix(i)

                M = D*M
                
                M_list.insert(0, D)
                D_list.insert(0, D)

            else:
                
                D = self.symInterfaceMatrix(i)
                P = self.symPropagationMatrix(i-1)

                M = P*D*M
                
                M_list.insert(0, P*D)
                D_list.insert(0, D)
                P_list.insert(0, P)

        r = M[1, 0]/M[0, 0]
        t = 1/M[0, 0]

        U = Matrix([t, 0])
        U = simplify(U)
        
        U_list = [U]

        for M_single in reversed(M_list):
            
            U = M_single*U
            U = simplify(U)

            U_list.insert(0, U)

        return M, M_list, D_list, P_list, U_list
    
    def createStringsFullMatrix(self, num_layers: int):
        
        num_interfaces = num_layers - 1
        
        r_str_list = []
        t_str_list = []
        phi_str_list = []
        
        for i in range(num_interfaces):
            
            interface_sub_str = '{sub0}{sub1}'.format(sub0 = i, sub1 = i+1)
            r_str = 'r' + interface_sub_str
            t_str = 't' + interface_sub_str

            r_str_list.append(r_str)
            t_str_list.append(t_str)

            if i != 0:
                
                layer_sub_str = '{sub0}'.format(sub0 = i)

                phi_str = 'phi' + layer_sub_str

                phi_str_list.append(phi_str)

        return r_str_list, t_str_list, phi_str_list 
    
    def fullMatrix(self, num_layers, r_list: list, t_list: list, phi_list: list):

        """ 
        Numerically evaluate the Transfer Matrix for n layers.
        
        Fresnel coefficients and phase factors for each interface/layer must be
        given.
        """

        M, M_list, D_list, P_list, U_list = self.symFullMatrix(num_layers)
        
        r_str_list, t_str_list, phi_str_list = self.createStringsFullMatrix(num_layers)

        r_zip = zip(r_str_list, r_list)
        t_zip = zip(t_str_list, t_list)
        phi_zip = zip(phi_str_list, phi_list)

        r_dict = dict(r_zip)
        t_dict = dict(t_zip)
        phi_dict = dict(phi_zip)

        full_dict = {**r_dict, **t_dict, **phi_dict}

        M = M.subs(full_dict)
        
        M_list = [np.array(simplify(i.subs(full_dict))).astype(np.complex128) for i in M_list]
        D_list = [np.array(simplify(i.subs(full_dict))).astype(np.complex128) for i in D_list]
        P_list = [np.array(simplify(i.subs(full_dict))).astype(np.complex128) for i in P_list]
        U_list = [np.array(simplify(i.subs(full_dict))).astype(np.complex128) for i in U_list]

        return M, D_list, P_list, U_list

if __name__ == "__main__":

    import numpy as np

    num_layers = 4

    TM = TransferMatrix()
    M, M_list, D_list, P_list, U_list = TM.symFullMatrix(num_layers)

    for i, U in enumerate(U_list):
    
        a = "a{i}".format(i = i) + " " + "(" + str(U_list[i][0]).replace(" ", "") + ")"
        n = "\n"
        b = "b{i}".format(i = i) + " " + "(" + str(U_list[i][1]).replace(" ", "") + ")"
    
        print(n)
        print(a)
        print(n)
        print(b)