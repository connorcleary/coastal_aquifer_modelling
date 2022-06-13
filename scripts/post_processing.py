import coastal_aquifer_model as cam
from pars import ModelParameters, load_parameters
import numpy as np

def find_mixing_volume(conc, pars, fraction=0.05):
    """
        Find the mixing zone volume

        Inputs:
            conc: 2D concentation array
            fraction: what fraction of salt water to consider 
            pars: ModelParameters object
        Outputs:
            mix:  volume of the mixing zone
    """
    mix = 0
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if 35.*fraction <= conc[i, j] <= 35*(1-fraction):
                mix += 1
    mixing_zone = pars.Ly*(pars.Lx/pars.ncol)*(pars.Lz/pars.nlay)*mix

    return mixing_zone


def find_toe_position(conc, pars, fraction=0.05):
    """
        Find the toe position

        Inputs:
            conc: 2D concentation array
            fraction: what fraction of salt water to consider 
            pars: ModelParameters object
        Outputs:
            toe: position of the toe
    """

    for i in range(pars.ncol-1, -1 , -1):       
        if conc[-1][i] > 35*fraction:
            return (i-pars.ncol*pars.offshore_proportion)*(pars.Lx/pars.ncol)


def find_mixing_centroid(conc, pars, fraction=0.05):
    """
        Find the centroid of the mixing zone

        Inputs:
            conc: 2D concentation array
            fraction: what fraction of salt water to consider 
            pars: ModelParameters object
        Outputs:
            centroid: centroid of the mixing
    """

    x_tot = 0 
    mix_n = 0

    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if 35.*fraction <= conc[i, j] <= 35*(1-fraction):
                mix_n += 1 
                x_tot += (j-pars.ncol*pars.offshore_proportion)*(pars.Lx/pars.ncol)

    centroid = x_tot/mix_n
    return centroid


def find_boundary_fluxes(conc, qx, qz, pars, fraction=0.05):
    """
        Find the boundary fluxes of the model

        Inputs: 
            conc: concetration array
            qx, qz: flux arrays
            pars: ModelParameters object
            fraction: fraction to consider saline
        Outputs:
            offshore_inflow_s: saline influx from the sea 
            offshore_inflow_f: fresh influx from the sea 
            offshore_outflow_s: saline outflux to the sea
            offshore_outflow_f: fresh outflux to the sea
            onshore_inflow_s: saline influx from the inland boundary
            onshore_inflow_f: fresh influx from the inland boundary
            onshore_outflow_s: saline outflux to the inland boundary
            onshore_outflow_f: fresh outflux to the inland boundary
    """

    delv = pars.Lz/pars.nlay

    offshore_inflow_s = 0 
    offshore_inflow_f = 0
    offshore_outflow_s = 0
    offshore_outflow_f = 0 
    onshore_inflow_s = 0
    onshore_inflow_f = 0
    onshore_outflow_s = 0
    onshore_outflow_f = 0

    # look at cells on the ends of the domain
    for k in range(pars.nlay):
        # offshore cells
        if k >= np.floor((pars.Lz-pars.sea_level)/delv): # if the cell is below sea level
            # if saline
            if conc[k, 0]>=35*fraction:
                # if influx
                if qx[k, 0]>0:
                    offshore_inflow_s += qx[k, 0]
                # if outflux
                else:
                    offshore_outflow_s -= qx[k, 0]
            # if fresh
            else:
                # if influx
                if qx[k, 0]>0:
                    offshore_inflow_f += qx[k, 0]
                # if outflux
                else:
                    offshore_outflow_f -= qx[k, 0]
    
        # onshore cells
        if k >= np.floor((pars.Lz-pars.sea_level-pars.h_b)/delv):
            # if saline
            if conc[k,-2]>=35*fraction:
                # if influx
                if qx[k,-2]<0:
                    onshore_inflow_s -= qx[k,-2]
                # if outflux
                else:
                    onshore_outflow_s += qx[k,-2]
            # if fresh
            else:
                # if influx
                if qx[k, -2]<0:
                    onshore_inflow_f -= qx[k,-2]
                # if outflux
                else:
                    onshore_outflow_f += qx[k,-2]
            
    # look at seafloor cells
    for i in range(int(pars.ncol*pars.offshore_proportion)):
        # if saline
        if conc[int((pars.Lz-pars.sea_level)/delv), i]>=35*fraction:
            # if influx
            if qz[int((pars.Lz-pars.sea_level)/delv), i]>0:
                offshore_inflow_s += qz[int((pars.Lz-pars.sea_level)/delv), i]
            # if outflux
            else:
                offshore_outflow_s -= qz[int((pars.Lz-pars.sea_level)/delv), i]
        # if fresh
        else:
            # if influx
            if qz[int((pars.Lz-pars.sea_level)/delv), i]>0:
                offshore_inflow_f += qz[int((pars.Lz-pars.sea_level)/delv), i]
            # if outflux
            else:
                offshore_outflow_f -= qz[int((pars.Lz-pars.sea_level)/delv), i]

    return offshore_inflow_s, offshore_inflow_f, offshore_outflow_s, offshore_outflow_f, \
           onshore_inflow_s, onshore_inflow_f, onshore_outflow_s, onshore_outflow_f


def find_mound(qx, pars):
    """
        Find the position of the mound, based on the smallest horizontal flux
            on the water table
        
        Inputs:
            qx: flux array
            pars: ModelParameters object
    """
    min_flux_value = np.inf
    min_flux_position = 0
    
    for k in range(pars.nlay):
        if float(np.max(qx[k,:])) != float(0):
            for i in range(pars.ncol):
                if float(qx[k,i]) != float(0) and np.abs(qx[k,i]) < min_flux_value:
                    min_flux_value = np.abs(qx[k,i])
                    min_flux_position = i
            break

    mound = (min_flux_position-pars.ncol*pars.offshore_proportion)*pars.Lx/pars.ncol

    return mound
                    