import numpy as np
import math
import random
from Analysis_useful_funcs import *

known_coords = {
'cartesian':{'name':['x', 'y', 'z'], 'metric':[['1','0','0'],['0','1','0'],['0','0','1']], 'sq_metric':[['1','0','0'],['0','1','0'],['0','0','1']]},
'cylindrical':{'name':['r', 'theta', 'z'], 'metric':[['1','0','0'],['0','r^2','0'],['0','0','1']], 'sq_metric':[['1','0','0'],['0','r','0'],['0','0','1']]},
'spherical':{'name':['r', 'theta', 'phi'], 'metric':[['1', '0', '0'],['0', 'r^2', '0'],['0', '0', 'r^2 sin^2(theta)']], 'sq_metric':[['1','0','0'],['0','r','0'],['0','0','r*sin(theta)']]}

}

def get_vol_term(h, idxs):
    '''inputs: h-vector of volume coefficients for non-orthonormal coordinate systems, idxs-which indexes from the volume coefficient vector we need
       outputs: vol_term-the volume term, made to look pretty (removing 1's and such
       This function returns the volume term that is necessary in getting the curl and divergence in non-orthonormal coordinate systems'''
       
    vol_term = ""
    
    for j in idxs:
        i = h[j]
        if is_float(i) and float(i) == 0:
            return "0"
        if not is_float(i) or float(i) != 1.0:
            vol_term += i
            if j +2 < len(h):
                vol_term += "*"
    return vol_term


def is_diagonal(Matrix):
    '''inputs: matrix-matrix which we want to check if is diagonal
       outputs: is_diag-boolean true or false
       This function just checks if a matrix is a diagonal matrix, WILL RETURN FALSE IF ANTIDIAGONAL'''
    
    matrix = np.asarray(Matrix)
    for i in range(len(matrix)):
        for j in range(len(matrix[i])):
            if i != j and float(matrix[i][j]) != 0.0:
                return False
                
                
    return True

def gradient(input_func, coord = 'Cartesian'):
    '''inputs: input_func-the function that I want the gradient of-should be a scalar with 1 dimension, coord-the RIGHT-HANDED coordinate system in which the input function is passed in (default to cartesian)
       outputs: grad-gradient of the input function
       This function takes the gradient of the input_func input function, assuming that it is in a right-handed coordinate system. The coordinate system can be passed in as well, defaults to the cartesian system'''
    
    tor = ""
    coord = coord.lower()
    
    if coord not in known_coords.keys():
        raise Exception("Gradient Function Unrecognized Coordinate System: "+str(coord))
    else:
        metric = np.asarray(known_coords[coord]['metric'])
        sq_metric = np.asarray(known_coords[coord]['sq_metric'])
        names = np.asarray(known_coords[coord]['name'])
        if not is_diagonal(metric):
            raise Exception("Non orthogonal Coordinate System: "+str(coord))
        else:
            diag = np.asarray(np.matrix(metric).diagonal())[0]
            
            h = np.asarray(np.matrix(sq_metric).diagonal())[0]
            
            tor += "[\n"
            for i in range(len(h)):
                tor += "(1/"+str(h[i])+")\td/d"+ str(names[i])+"["+str(input_func)+"]    e_{"+str(names[i])+"}"
                if i + 1 < len(h):
                    tor += "  ,\n"
            tor += "\n]"
            
    return tor

    
    

def divergence(input_func, coord = 'Cartesian'):
    '''inputs: input_func-the function that I want the divergence of-should be a vector with 3 dimensions, coord-the RIGHT-HANDED coordinate system in which the input function is passed in (default to cartesian)
       outputs: grad-gradient of the input function
       This function takes the divergence of the input_func input function, assuming that it is in a right-handed coordinate system. The coordinate system can be passed in as well, defaults to the cartesian system'''
    
    tor = ""
    coord = coord.lower()
    if len(input_func) != 3:
        raise Exception("Input function to gradient is not 3 dimensional")
    if coord not in known_coords.keys():
        raise Exception("Gradient Function Unrecognized Coordinate System: "+str(coord))
    else:
        metric = np.asarray(known_coords[coord]['metric'])
        sq_metric = np.asarray(known_coords[coord]['sq_metric'])
        names = np.asarray(known_coords[coord]['name'])
        if not is_diagonal(metric):
            raise Exception("Non orthogonal Coordinate System: "+str(coord))
        else:
            diag = np.asarray(np.matrix(metric).diagonal())[0]
            h = np.asarray(np.matrix(sq_metric).diagonal())[0]
            
            
            idxs = [0,1,2]
            vol_term = get_vol_term(h = h, idxs = idxs)
            for i in range(len(h)):
                if len(vol_term) > 0:
                    tor += "(1/"+vol_term+") * d/d"+ str(names[i])+"["+vol_term+"*"+str(input_func[i])+"]"
                else:
                    tor += "d/d"+ str(names[i])+"["+str(input_func[i])+"]"
                if i + 1 < len(h):
                    tor += "  +\n"
            
            
    return tor
    


def curl(input_func, coord = 'Cartesian'):
    '''inputs: input_func-the function that I want the curl of-should be a vector with 3 dimensions, coord-the RIGHT-HANDED coordinate system in which the input function is passed in (default to cartesian)
       outputs: grad-gradient of the input function
       This function takes the divergence of the input_func input function, assuming that it is in a right-handed coordinate system. The coordinate system can be passed in as well, defaults to the cartesian system'''
    
    tor = ""
    coord = coord.lower()
    if len(input_func) != 3:
        raise Exception("Input function to gradient is not 3 dimensional")
    if coord not in known_coords.keys():
        raise Exception("Gradient Function Unrecognized Coordinate System: "+str(coord))
    else:
        metric = np.asarray(known_coords[coord]['metric'])
        sq_metric = np.asarray(known_coords[coord]['sq_metric'])
        names = np.asarray(known_coords[coord]['name'])
        if not is_diagonal(metric):
            raise Exception("Non orthogonal Coordinate System: "+str(coord))
        else:
            diag = np.asarray(np.matrix(metric).diagonal())[0]
            h = np.asarray(np.matrix(sq_metric).diagonal())[0]
            
            #(curl F)1 = (1/h2*h3)[d(h3F3)/du2 - d(h2F2)/du3]
            #(curl F)2 = (1/h1*h3)[d(h1F1)/du3 - d(h3F3)/du1]
            #(curl F)3 = (1/h1*h2)[d(h2F2)/du1 - d(h1F1)/du2]
            
            #doing (curl F)1
            idxs = [1,2]
            vol_term = get_vol_term(h = h, idxs = idxs)
            if len(vol_term) > 0:
                curl_F_1 = "(1/"+str(vol_term)+")\t[d("+str(h[2])+"*"+str(input_func[2])+")/d("+str(names[1])+") - d("+str(h[1])+"*"+str(input_func[1])+")/d("+str(names[2])+")]"
            else:
                curl_F_1 = "[d("+str(h[2])+"*"+str(input_func[2])+")/d("+str(names[1])+") - d("+str(h[1])+"*"+str(input_func[1])+")/d("+str(names[2])+")]"
            
            
            idxs = [0,2]
            vol_term = get_vol_term(h = h, idxs = idxs)
            if len(vol_term) > 0:
                curl_F_2 = "(1/"+str(vol_term)+")\t[d("+str(h[0])+"*"+str(input_func[0])+")/d("+str(names[2])+") - d("+str(h[2])+"*"+str(input_func[2])+")/d("+str(names[0])+")]"
            else:
                curl_F_2 = "[d("+str(h[0])+"*"+str(input_func[0])+")/d("+str(names[2])+") - d("+str(h[2])+"*"+str(input_func[2])+")/d("+str(names[0])+")]"
               
               
            idxs = [0,1]
            vol_term = get_vol_term(h = h, idxs = idxs)
            if len(vol_term)>0:
                curl_F_3 = "(1/"+str(vol_term)+")\t[d("+str(h[1])+"*"+str(input_func[1])+")/d("+str(names[0])+") - d("+str(h[0])+"*"+str(input_func[0])+")/d("+str(names[1])+")]"
            else:
                curl_F_3 = "[d("+str(h[1])+"*"+str(input_func[1])+")/d("+str(names[0])+") - d("+str(h[0])+"*"+str(input_func[0])+")/d("+str(names[1])+")]"

            tor = "[\n"+curl_F_1+"  e_{"+str(names[0])+"}   ,\n"+curl_F_2+"  e_{"+str(names[1])+"}   ,\n"+curl_F_3+"  e_{"+str(names[2])+"}\n]"
            
            
    return tor
