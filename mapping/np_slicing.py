'''
Created on Mar 31, 2015

@author: ayan
'''
import numpy as np


if __name__ == '__main__':
    
    t1 = (10, 12, 13, 14)
    t2 = (20, 22, 23, 15)
    t3 = (41, 68, 32, 19)
    x = np.array([t1, t2, t3])
    s = np.s_[:, 1:3]
    print(s)
    print(isinstance(s, tuple))
    print(x)[s]