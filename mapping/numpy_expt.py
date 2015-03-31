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
    
    print('\n')
    
    s1 = (1, 2)
    s2 = (3, 4)
    s3 = (5, 6)
    b = (s1, s2, s3)
    y = np.array(b)
    print(y)
    yt = np.transpose(y)
    yt_trim_high = yt[:, 1:]
    yt_trim_low = yt[:, :-1]
    yt_avg = 0.5 * (yt_trim_high + yt_trim_low)
    y_avg = np.transpose(yt_avg)
    print(y_avg)
    
    print('\n')
    
    
    