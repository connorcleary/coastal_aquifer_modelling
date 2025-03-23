import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters

def _create_temporal_discretization(pars, sea_level_rise, rise_length, rise_type, rise_rate, initial_conditions):
    perlen = []
    nstep = []
    steady = []
    if initial_conditions is None:
        perlen.append(pars.perlen)
        nstep.append(pars.perlen/pars.dt)
        steady.append(True)
        if sea_level_rise == 0:
            return perlen, nstep, steady
    if rise_type is 'step':
        perlen.append(rise_length)
        nstep.append(np.round(rise_length/pars.dt))
        steady.append(False)
        return perlen, nstep, steady
    if rise_type is 'linear':
        for i in range(int(rise_length/pars.dt)):
            perlen.append(pars.dt)
            nstep.append(1)
            steady.append(False)
        return perlen, nstep, steady        
    
def _create_cell_groups(pars, delr, delv, delc):
    # define cell groups
    inactive_cells = []
    offshore_boundary_cells = []
    onshore_boundary_cells = []
    #  surface_boundary_cells = []
    wetland_cells = []

    # add inactive cells
    for i in range(int(pars.ncol * pars.offshore_proportion)):
        for j in range(pars.nrow):
            for k in range(0, int((pars.Lz - pars.sea_level) / delv)):
                inactive_cells.append([k, j, i])

    # add cells on ends of domain
    for k in range(pars.nlay):
        for j in range(pars.nrow):
            if k >= np.floor((pars.Lz - pars.sea_level) / delv):  # if the cell is below sea level
                offshore_boundary_cells.append([k, j, 0])

            if k >= np.floor((pars.Lz - pars.sea_level - pars.h_b) / delv):
                onshore_boundary_cells.append([k, j, pars.ncol - 1])

    # add the seafloor
    for i in range(int(pars.ncol * pars.offshore_proportion)):
        for j in range(pars.nrow):
            offshore_boundary_cells.append([int((pars.Lz - pars.sea_level) / delv), j, i])

    # add wetland cells
    if pars.x_w != 0:
        for j in range(pars.nrow):
            for i in range(int((pars.offshore_proportion * pars.Lx + pars.x_w) / delr),
                           int((pars.offshore_proportion * pars.Lx + pars.x_w + pars.Lx_w) / delr)):
                for k in range(int((pars.Lz - pars.sea_level - pars.h_w) / delv) + 1,
                               int((pars.Lz - pars.sea_level - pars.z_w) / delv)):
                    wetland_cells.append([k, j, i])
                for k in range(0, int((pars.Lz - pars.sea_level - pars.h_w) / delv) + 1):
                    inactive_cells.append([k, j, i])

    # create ibound array
    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)
    for cell in inactive_cells:
        ibound[cell[0], cell[1], cell[2]] = 0
    else:
        for cell in onshore_boundary_cells + offshore_boundary_cells:
            ibound[cell[0], cell[1], cell[2]] = -1
            
    return inactive_cells, offshore_boundary_cells, onshore_boundary_cells, wetland_cells, ibound

def _create_stress_period_data(pars, nper, sea_level_rise, rise_length, rise_type, rise_rate, initial_conditions,
                               wetland_cells, onshore_boundary_cells, offshore_boundary_cells):
    drn_spd = {per: [] for per in range(nper)}
    chd_spd = {per: [] for per in range(nper)}
    ssm_spd = {per: [] for per in range(nper)}
    oc_spd = {}
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    for kstp in range(0, int(pars.perlen / pars.dt), pars.frequency):
        oc_spd[(0, kstp)] = ["save head", "save budget"]

    for cell in wetland_cells:
        for per in range(nper):
            drn_spd[per].append([cell[0], cell[1], cell[2], pars.sea_level+pars.h_w, pars.drain_conductance])
            ssm_spd[per].append([cell[0], cell[1], cell[2], 0.0, itype["DRN"]])

    for cell in onshore_boundary_cells:
        for per in range(nper):
            ssm_spd[per].append([cell[0], cell[1], cell[2], 0, itype["BAS6"]])
            chd_spd[per].append([cell[0], cell[1], cell[2], pars.sea_level+pars.h_b, pars.sea_level+pars.h_b])

    if initial_conditions is None:
        for cell in offshore_boundary_cells:
            ssm_spd[0].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
            chd_spd[0].append([cell[0], cell[1], cell[2], pars.sea_level, pars.sea_level])
            for kstp in range(0, int(pars.perlen / pars.dt), pars.frequency):
                oc_spd[(0, kstp)] = ["save head", "save budget"]

    if rise_type == 'step':
        for cell in offshore_boundary_cells: #+ _get_extra_chd_cells(pars, sea_level_rise):
            ssm_spd[1].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
            chd_spd[1].append([cell[0], cell[1], cell[2], pars.sea_level+sea_level_rise, pars.sea_level+sea_level_rise])
            for kstp in range(0, int(rise_length / pars.dt), pars.frequency):
                oc_spd[(1, kstp)] = ["save head", "save budget"]

    if rise_type == 'linear':
        for per in range(1, nper):
            sea_level_rise_increment = per/nper*sea_level_rise
            for cell in offshore_boundary_cells: #+ _get_extra_chd_cells(pars, sea_level_rise_increment):
                ssm_spd[per].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
                chd_spd[per].append([cell[0], cell[1], cell[2], pars.sea_level + sea_level_rise_increment,
                                   pars.sea_level + sea_level_rise_increment])
            if per%pars.frequency == 0 or per == nper-1:
                oc_spd[(per, 0)] = ["save head", "save budget"]

    return drn_spd, chd_spd, oc_spd, ssm_spd

def build_model(pars, sea_level_rise=0, rise_type='linear', rise_length=100*365, rise_time_series=None, initial_conditions=None):
    '''
        A function to build a coastal aquifer model.
        :param pars: parameters defining the geometry of the model
        :param steady: bool, whether the model is steady or not
        :param sea_level_rise: float, int, how much sea level rise to apply
        :param rise_length: float, int, time to apply sea-level rise over (days)
        :param rise_type: string, how the sea level rise is applied, one 'step', 'linear' # todo implement for custom
        :param rise_time_series: if rise_type is 'custom', then this a time series showing the curve
        :param initial_conditions: string, name of the model to get initial conditions from.
        :return:
    '''

    # create model workspace
    model_ws = f".\\model_files\\{pars.name}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    # create base seawat model
    swt = flopy.seawat.Seawat(pars.name, model_ws=model_ws, exe_name=pars.exe_path)

    # calculate cell dimension
    delr = pars.Lx/pars.ncol
    delc = pars.Ly/pars.nrow
    delv = pars.Lz/pars.nlay

    # define top and bottoms
    assert pars.Lz-pars.sea_level >= sea_level_rise, "Model top must be higher than sea-level"
    top = pars.Lz
    botm = np.linspace(top-delv, 0, pars.nlay)

    # something I've copied
    ipakcb = 53

    perlen, nstep, steady = _create_temporal_discretization(pars, sea_level_rise, rise_length, rise_type, rise_time_series, initial_conditions,
                                                           )

    # define discretization package
    dis = flopy.modflow.ModflowDis(
            model=swt,
            nlay=pars.nlay,
            nrow=pars.nrow,
            ncol=pars.ncol,
            nper=len(perlen),
            itmuni=4, # four for days
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=perlen,
            nstp=nstep,
            steady=steady
        )

    inactive_cells, offshore_boundary_cells, onshore_boundary_cells, wetland_cells, ibound = _create_cell_groups(pars, delr, delv, delc)

    # define starting heads
    if initial_conditions is None:
        strt = pars.Lz*np.ones((pars.nlay, pars.nrow, pars.ncol))
    else:
        raise NotImplementedError

    # create basic package
    bas = flopy.modflow.ModflowBas(
            model=swt, 
            ibound=ibound, 
            strt=strt
        )

    # define layer types and wetting, this one I'm not sure about
    laytyp=np.ones(pars.nlay)
    laytyp[0] = 1
    laywet=np.ones(pars.nlay)
    laywet[0] = 1

    # create layer property flow package
    lpf = flopy.modflow.ModflowLpf(
            swt, 
            hk=pars.K, 
            vka=pars.anis, 
            ipakcb=ipakcb, 
            laytyp=laytyp, 
            laywet=laywet,
            ss=pars.ss, # not sure about these ones
            sy=pars.sy,
            layvka=1,
        )

    # create solver package
    pcg = flopy.modflow.ModflowPcg(
            swt, 
            hclose=1.0e-5, 
            npcond=1, 
            mxiter=500
        )

    drn_spd, chd_spd, oc_spd, ssm_spd = _create_stress_period_data(pars, len(perlen), sea_level_rise, rise_length, rise_type, rise_time_series, initial_conditions,
                                                                   wetland_cells, onshore_boundary_cells, offshore_boundary_cells)

    if len(wetland_cells) > 0:
        drn = flopy.modflow.ModflowDrn(
            model=swt,
            stress_period_data=drn_spd,
            ipakcb=ipakcb
        )

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # define constant head package
    chd = flopy.modflow.ModflowChd(
            model=swt, 
            stress_period_data=chd_spd,
            ipakcb=ipakcb
        )

    # create recharge package
    rch = flopy.modflow.ModflowRch(
            model=swt,
            rech=pars.W_net, #rech,
            ipakcb=ipakcb
        )

    # set starting concentrations
    if initial_conditions is None:
        sconc = 0.0*np.ones((pars.nlay, pars.nrow, pars.ncol))
        sconc[:, :, 0] = 35.0
    else:
        raise NotImplementedError

    # define basic transport package
    btn = flopy.mt3d.Mt3dBtn(
            swt,
            nprs=-pars.frequency,
            prsity=pars.n,
            sconc=sconc,
            chkmas=False,
            nprobs=10,
            nprmas=10,
            dt0=pars.dt
        )

    # define advection package
    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)

    # define dispersion package
    dsp = flopy.mt3d.Mt3dDsp(
            swt, 
            al=pars.alpha_L, 
            trpt=pars.alpha_anisT, 
            trpv=pars.alpha_anisV, 
            dmcoef=pars.diff,
        )

    # define transport solver
    gcg = flopy.mt3d.Mt3dGcg(
            model=swt, 
            iter1=500, 
            mxiter=1, 
            isolve=2, 
            cclose=1e-5
        )

    # find number of sinks and sources
    mxss = int(np.ceil(2*pars.nlay*pars.nrow + 
                        pars.nrow*pars.ncol*pars.offshore_proportion+1 +
                        pars.nrow*pars.ncol+pars.Lx_w*pars.nrow*(pars.x_w!=0)))

    # define source sink mixing package
    ssm = flopy.mt3d.Mt3dSsm(
            model=swt,
            stress_period_data=ssm_spd,
            mxss=mxss
        )

    vdf = flopy.seawat.SeawatVdf(
            swt,
            iwtable=0,
            densemin=0,
            densemax=0,
            denseref=1000.0,
            denseslp=0.7143,
            firstdt=pars.dt,
        )
    # write input
    swt.write_input() 

    return swt


def run_model(swt):
    """
        A function to run the seawat model

        Inputs: 
            swt: model object
        Outputs:
            None
    """
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")


def extract_results(name):
    """
        Open model results from binary files

        Inputs:
            name: name of model/realization/scenario
        Outputs:
            head: head matrix [nstp, nlay, nrow, ncol]
            qx: longitudinal flux matrix [nstp, nlay, nrow, ncol]
            qy: transverse flux matrix matrix [nstp, nlay, nrow, ncol]
            qz: vertical flux matrix matrix [nstp, nlay, nrow, ncol]
            concentration: concentration matrix [nstp, nlay, nrow, ncol]
    """
    pars = load_parameters(name)
    name = pars.name
    model_ws = f".\\model_files\\{name}"
    nstp = pars.perlen/pars.dt

    # open binary files
    ucnobj = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws, f'{name}.hds'))

    # get head and concentration data
    concentration = ucnobj.get_alldata()[:]
    head = headobj.get_alldata()[:]
    
    # select every n items
    times = ucnobj.get_times()
    concentration = concentration

    qx = np.zeros_like(concentration)
    qy = np.zeros_like(concentration)
    qz = np.zeros_like(concentration)

    # get fluxes
    for t in range(qx.shape[0]):
        qx[t] = cbbobj.get_data(text="flow right face", totim=times[t])[0]
        if pars.nrow > 1:
            qy[t] = cbbobj.get_data(text="flow front face", totim=times[t])[0]
        qz[t] = cbbobj.get_data(text="flow lower face", totim=times[t])[0]

    save_results(name, concentration, head, qx, qy, qz)
    return concentration, head, qx, qy, qz


def save_results(name, concentration, head, qx, qy, qz):
    """
        Save extracted results to a .npy file

        Inputs:
            name: model name
            concentration, head etc. : numpy arrays of model outputs
        Outputs:
            None
    """
    ws = os.path.join(f'.\\results\\{name}')
    if not os.path.exists(ws):
        os.makedirs(ws)

    with open(os.path.join(ws, f"qx.npy"), 'wb') as f: np.save(f, np.array(qx))
    with open(os.path.join(ws, f"qy.npy"), 'wb') as f: np.save(f, np.array(qy))
    with open(os.path.join(ws, f"qz.npy"), 'wb') as f: np.save(f, np.array(qz))
    with open(os.path.join(ws, f"head.npy"), 'wb') as f: np.save(f, np.array(head))
    with open(os.path.join(ws, f"concentration.npy"), 'wb') as f: np.save(f, np.array(concentration))


def load_results(name):
    """
        Load extracted results from .npy files

        Inputs:
            name: name of the model
        Outputs:
            concentration, head... : numpy matrices of results
    """
    ws = os.path.join(f'.\\results\\{name}')

    with open(os.path.join(ws, f"qx.npy"), 'rb') as f: qx = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qy.npy"), 'rb') as f: qy = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qz.npy"), 'rb') as f: qz = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"head.npy"), 'rb') as f: head = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"concentration.npy"), 'rb') as f: concentration = np.load(f, allow_pickle=True)

    return concentration, head, qx, qy, qz, 