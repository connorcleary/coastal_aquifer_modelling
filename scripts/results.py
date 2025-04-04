from pars import ModelParameters, load_parameters
import numpy as np
import coastal_aquifer_model as cam
import matplotlib.pyplot as plt
import plot_helpers as plth
import os
import post_processing as proc

def plot_salinity(pars, ax, conc, qx, qz, cmap="viridis", row=0, x_step=10, z_step=5, width=0.002, arrow_c="white"):

    x = np.linspace(-pars.Lx*pars.offshore_proportion, pars.Lx-pars.Lx*pars.offshore_proportion, pars.ncol)
    y = np.linspace(-pars.sea_level, pars.Lz-pars.sea_level, pars.nlay)

    concentration_array = conc[:, row, :]
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if concentration_array[i, j] == np.float32(1.e30):
                concentration_array[i, j] = np.nan

    conccm = ax.pcolormesh(x, y, np.flipud(concentration_array),
                               cmap=cmap, vmax=35, vmin=0)

    # plot arrows
    X, Y = np.meshgrid(x[::x_step], y[::z_step])
    ax.quiver(X, Y, np.flipud(qx[::z_step, row, ::x_step]), np.flipud(qz[::z_step, row, ::x_step]),
                  color=arrow_c, width=width)
    ax.set_box_aspect(0.25)
    return ax



def plot_results(name, timestep=-1, row=0, return_axs=False, figsize=(12,6),
    cmap="viridis", arrow_c="white", aspect=8, x_step=10, z_step=10, width=0.002, fmt="%3.2f"):
    """
        Plot the head, concentrations and fluxes

        Inputs:
            name: model to plot
            timestep: timestep to plot, default is -1 (the last one)
            row: row to plot
            return_axs: flag whether to return the axes objects
            figsize: figure dimensions (inches)
            cmap: colormap
            arrow_c: color of arrows
            aspect: vertical exageration 
            vector_T: spaces between vector arrows
            width: arrow width
            fmt: format of contour labels
        Outputs:
            axs: axes objects (optional)
    """

    # load parameters and results
    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)

    f, axs = plt.subplots(2, 1, figsize=figsize)

    # set up x and y arrays, to be distance above sea level and distance onshore
    x = np.linspace(-pars.Lx*pars.offshore_proportion, pars.Lx-pars.Lx*pars.offshore_proportion, pars.ncol)
    y = np.linspace(-pars.sea_level, pars.Lz-pars.sea_level, pars.nlay)

    # select relevent slice in time and the alongshore direction, and set values above the water table as nan
    concentration_array = concentration[timestep,:,row,:]
    head_array = head[timestep,:,row,:] - pars.sea_level*np.ones_like(head[timestep,:,row,:])
    for i in range(pars.nlay):
        for j in range(pars.ncol):
            if concentration_array[i, j] == np.float32(1.e30):
                head_array[i,j] = np.nan
                concentration_array[i, j] = np.nan

    # plot head colormesh
    headcm = axs[0].pcolormesh(x, y, np.flipud(head_array), 
            	                cmap=cmap, vmax=np.nanmax(head_array[:, 1:]), vmin=np.min([0, pars.h_b]))

    # plot head contours
    hc = axs[0].contour(x, y, np.flipud(head_array), colors=arrow_c, 
                        levels=np.linspace(np.min([0, pars.h_b]), np.nanmax(head_array[:, 1:]), 15))

    # label contours
    axs[0].clabel(hc, hc.levels, inline=True, fontsize=10, fmt=fmt)

    # plot concentration colormesh
    conccm = axs[1].pcolormesh(x, y, np.flipud(concentration_array), 
            	                cmap=cmap, vmax=35, vmin=0)

    # plot arrows
    X, Y = np.meshgrid(x[::x_step],y[::z_step])
    axs[1].quiver(X, Y,
                    np.flipud(qx[timestep,::z_step,row,::x_step]),
                    np.flipud(qz[timestep,::z_step,row,::x_step]),
                    color=arrow_c, width=width)

    axs[0].set_aspect(aspect)
    axs[1].set_aspect(aspect)
    axs[0].set_title("Head")
    axs[1].set_title("Salinity")
    axs[0].set_ylabel("Height above sealevel (m)")
    axs[1].set_ylabel("Height above sealevel (m)")
    axs[1].set_xlabel("Distance onshore (m)")

    f.suptitle(f"Head and salinity distributions for {name}")

    headcb = plt.colorbar(headcm, shrink=1, ax=axs[0])
    conccb = plt.colorbar(conccm, shrink=1, ax=axs[1])
    headcb.ax.set_title('Head (m)', fontsize = 'small')
    conccb.ax.set_title('Salinity (kg/m^3)', fontsize = 'small')
    

    # return axs objects if necessary
    if return_axs:
        return axs
    else:
        ws = os.path.join(f'.\\figures\\{name}')
        if not os.path.exists(ws):
            os.makedirs(ws)
        plt.savefig(f"{ws}\\head_and_concentration", dpi=300)


def plot_evolutions(name, row=0, fraction = 0.01, return_axs=False, figsize=(18,6), interval=20):
    """
        Plot evolutions of metrics, to check a steady state has 
        been reached

        Inputs:
            name: name of the model
            row: row to plot
            fraction: fraction of saltwater to consider as the mixing 
                zone i.e. 35*fraction <= mixing zone <= 35*(1-fraction)
            return_axs: optionally return the axis
            figsize: size of figure in inches
            interval: number of timesteps between metric calculations
    """
    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)
    nstp = int(pars.perlen/pars.dt)
    # set as 1% saltwater

    # create arrays
    toe = np.zeros(int(nstp/interval))
    mixing_volume = np.zeros(int(nstp/interval))
    centroid = np.zeros(int(nstp/interval))

    # select every 20th step
    times = np.linspace(0, nstp-1, int(nstp/interval), dtype=int)

    for tdx, t in enumerate(times):
        toe[tdx] = proc.find_toe_position(concentration[t, :, row, :], pars, fraction)
        mixing_volume[tdx] = proc.find_mixing_volume(concentration[t, :, row, :], pars, fraction)
        centroid[tdx] = proc.find_mixing_centroid(concentration[t, :, row, :], pars, fraction)


    f, axs = plt.subplots(3, 1, sharex=True, figsize=figsize)

    # plot lines
    axs[0].plot(pars.dt*times, toe)
    axs[1].plot(pars.dt*times, mixing_volume)
    axs[2].plot(pars.dt*times, centroid)

    axs[0].set_ylabel("distance onshore (m)")
    axs[1].set_ylabel("volume (m^3)")
    axs[2].set_ylabel("distance onshore (m)")

    axs[2].set_xlabel("time (days)")

    axs[0].set_title("Toe position")
    axs[1].set_title("Mixing zone volume")  
    axs[2].set_title("Mixing zone centroid")

    f.suptitle(f"Evolution of metrics for {name}")

    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)
        
    plt.savefig(f"{ws}\\metric_evolutions", dpi=300)

    if return_axs: return axs


def plot_boundary_concentration(name, return_axs=False, figsize=(6,6), row = 0):
    """
        Plot the concentrations in the cells along the right hand boundary

        Inputs:
            name: model name
            return_axs: flag whether to return axes
            figsize: figure size in inches
            row: model row
    """

    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)

    # find mound position
    mound = proc.find_mound(qx[-1,:,row,:], pars)
    mound_col = np.int((mound/pars.Lx+pars.offshore_proportion)*pars.ncol)

    # truncate to just include layers below the inland boundary
    concentration_mound = np.copy(concentration[-1,int(pars.nlay/pars.Lz*(pars.Lz-pars.sea_level-pars.h_b)):,row,mound_col])
    concentration_edge = np.copy(concentration[-1,int(pars.nlay/pars.Lz*(pars.Lz-pars.sea_level-pars.h_b)):,row,-1])

    y=np.linspace(-pars.sea_level, pars.h_b, len(concentration_edge))

    inland_lay = int(pars.nlay/pars.Lz*(pars.Lz-pars.sea_level-pars.h_b))+1
    # find isoclors at the mound column
    concentration_mound[-1] = 35
    lay_isoclor1 = np.atleast_1d(np.argmax(concentration_mound>0.35))[0] + inland_lay
    lay_isoclor5 = np.atleast_1d(np.argmax(concentration_mound>1.75))[0] + inland_lay
    lay_isoclor10 = np.atleast_1d(np.argmax(concentration_mound>3.5))[0] + inland_lay

    f, ax = plt.subplots(figsize=figsize)
    ax.plot(np.flipud(concentration_edge), y)
    ax.set_ylabel("Distance above sea level (m)")
    ax.set_xlabel("Salinity (kg/m^3)")
    ax.set_title(f"Salinity at inland boundary for {name}")

    # plot isoclors
    delv = pars.Lz/pars.nlay
    ax.axhline(pars.Lz-pars.sea_level-(lay_isoclor1)*delv, c='b', alpha=0.5, zorder = -1, linestyle=':', label=r"mound 1% isoclor")
    ax.axhline(pars.Lz-pars.sea_level-(lay_isoclor5)*delv, c='g', alpha=0.5, zorder = -1, linestyle='-', label=r"mound 5% isoclor")
    ax.axhline(pars.Lz-pars.sea_level-(lay_isoclor10)*delv, c='r', alpha=0.5, zorder = -1, linestyle='--', label=r"mound 10% isoclor")

    ax.set_xlim([-5, 40])
    ax.set_ylim([-pars.sea_level-0.5, pars.Lz-pars.sea_level+0.5])
    ax.legend()

    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)
        
    plt.savefig(f"{ws}\inland_boundary_salinity", dpi=300)

    if return_axs: return ax

def plot_sea_level_rise_results(name):
    f, axs = plt.subplots(4, 3, figsize=(10,6))
    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)
    axs = axs.flatten()
    if pars.initial_conditions is None:
        i_init = int(pars.perlen/pars.dt/pars.frequency)
        plot_salinity(pars, axs[0], concentration[i_init], qx[i_init], qz[i_init])
        axs[0].set_title("Steady state")
    else:
        raise NotImplementedError("still need to set up using different initial conditions")

    for i in range(1, 11):
        step = int(pars.rise_length/10/(pars.dt*pars.frequency))
        plot_salinity(pars, axs[i], concentration[i_init+i*step], qx[i_init+i*step], qz[i_init+i*step])
        axs[i].set_title(f"{i*pars.rise_length/10/365:1.0f} years")

    axs[-1].set_axis_off()

    plt.tight_layout()
    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    plt.savefig(f"{ws}\concentration_snapshots.png", dpi=600)

def plot_drain_discharge_and_salinity(name):

    f, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    pars = load_parameters(name)
    if pars.initial_conditions is not None:
        raise NotImplementedError

    results = cam.get_drain_results(name)
    idx = results['times'] > results['times'][-1] - pars.rise_length
    times = (results['times'][idx] - (results['times'][-1]- pars.rise_length))/365
    discharges = results['discharge'][idx]
    concentrations = results['conc'][idx]

    axs[0].plot(times, -discharges)
    axs[0].set_ylabel("Discharge [m^2/day]")
    axs[1].plot(times, concentrations)
    axs[1].set_ylabel("Concentration [PSU]")
    axs[1].set_xlabel("Time [years]")

    ws = os.path.join(f'.\\figures\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    plt.savefig(f"{ws}\drain_flow_and_concentration.png", dpi=600)



def save_metrics(name, row=0, fraction=0.05):
    """
        Find and save metrics for the final timestep. Saves to text file

        Inputs:
            name: name of model
            row: row of model
            fraction: fraction to consider saline
        Outputs: 
            none
    """

    pars = load_parameters(name)
    concentration, head, qx, qy, qz = cam.load_results(name)

    offshore_inflow_s, offshore_inflow_f, offshore_outflow_s, offshore_outflow_f, \
    onshore_inflow_s, onshore_inflow_f, onshore_outflow_s, onshore_outflow_f = \
    proc.find_boundary_fluxes(concentration[-1,:,row,:], \
        qx[-1,:,row,:], qz[-1,:,row,:], pars, fraction=fraction)

    toe = proc.find_toe_position(concentration[-1, :, row, :], pars, fraction)
    mixing_volume = proc.find_mixing_volume(concentration[-1, :, row, :], pars, fraction)
    centroid = proc.find_mixing_centroid(concentration[-1, :, row, :], pars, fraction)

    mound = proc.find_mound(qx[-1,:,row,:], pars)

    metrics = [toe, centroid, mixing_volume, mound, offshore_inflow_s, 
               offshore_inflow_f, offshore_outflow_s, offshore_outflow_f,
               onshore_inflow_s, onshore_inflow_f, onshore_outflow_s, 
               onshore_outflow_f]

    strings = ["Position of toe (distance onshore)", 
        	   "Centroid of mixing zone (distance onshore)",
               "Area of mixing zone",
               "Mound position (distance onshore)",
               "Saline inflow from the sea",
               "Fresh inflow from the sea",
               "Saline outflow to the sea",
               "Fresh outflow to the sea",
               "Saline inflow from the inland boundary",
               "Fresh inflow from the inland boundary",
               "Saline outflow to the inland boundary",
               "Fresh outflow to the inland boundary"]

    units = ["m", "m", "m^3", "m", "m^3/s","m^3/s","m^3/s",
             "m^3/s","m^3/s","m^3/s","m^3/s","m^3/s"]

    ws = os.path.join(f'.\\results\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    with open(f"{ws}\\metrics_{fraction}.txt", "w") as f:

        f.write(f"Metrics for {name} with fraction={fraction} \n\n")
        for (metric, string, unit) in zip(metrics, strings, units):

            f.write(f"{string}: {metric} {unit}\n")

        if pars.x_w != 0:
            f.write(f"Concentration flux to wetland: {proc.find_wetland_flux(concentration[-1,:,row,:], qx[-1,:,row,:], qz[-1,:,row,:], pars)} kg m^-1 day^-1")

    concentration_b = concentration[-1, :, row, -1]
    # cell centres
    depths = np.linspace(pars.Lz-pars.sea_level-(pars.Lz/(2*pars.nlay)), -pars.sea_level+(pars.Lz/(2*pars.nlay)), pars.nlay)
    lays = np.linspace(0,109, 110)

    headers = ["layer", "masl", "C"]
    results = np.vstack((headers, np.stack((lays, depths, concentration_b)).T))
    np.savetxt(f"{ws}\\boundary_concentrations.csv", results, delimiter=",", fmt ='% s')