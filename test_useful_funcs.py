import numpy as np
import math
import useful_funcs
from useful_funcs import *




to_clean = [0.0001, 0.001, 100.00001, 10000000.01]
for t in to_clean:
    print(str_clean_number(t, num_digits = 6))
    print("")
    
