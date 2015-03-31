'''
Created on Mar 31, 2015

@author: ayan
'''
import numpy as np


if __name__ == '__main__':
    
    t1 = (1, 2, 3)
    t2 = (10, 20, 30)
    a = [t1, t2]
    x = np.array(a)
    x_trim_high = x[:, 1:]
    x_trim_low = x[:, :-1]
    print(x_trim_high)
    print(x_trim_low)
    x_row_avg = 0.5 * (x_trim_high + x_trim_low)
    print(x_row_avg)