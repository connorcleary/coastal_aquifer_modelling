import matplotlib.pyplot as plt
import numpy as np

def sample_grid_data(array, vstp, hstp):

    if len(array.shape) == 2:
        rows, cols = array.shape
    else:
        cols = array.shape[0]
        rows = 1

    if len(array.shape) == 2:
        new_array = np.zeros((((rows-1)//vstp+1), (cols-1)//hstp+1))
    else:
        new_array = np.zeros((cols-1)//hstp+1)
        
    for n in range(rows):
        for m in range(cols):  
            if n%vstp == 0 and m%hstp == 0:
                if len(array.shape) == 2:
                    new_array[n//vstp, m//hstp] = array[n,m]
                else:
                    new_array[m//hstp] = array[m]

    return new_array