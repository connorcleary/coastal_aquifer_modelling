import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
from pars import ModelParameters, load_parameters


def build_steady_model(pars):
    '''
        A function to build a coastal aquifer model.

        Input: 
            pars: parameters object

        Outputs:
            swt: groundwater flow and transport model
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
    top = pars.Lz
    botm = np.linspace(top-delv, 0, pars.nlay)

    # something I've copied
    ipakcb = 53

    # define discretization package
    dis = flopy.modflow.ModflowDis(
            model=swt,
            nlay=pars.nlay,
            nrow=pars.nrow,
            ncol=pars.ncol,
            nper=1,
            itmuni=4, # four for days
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=pars.perlen,
            nstp=pars.perlen/pars.dt,
            steady=True
        )

    # define cell groups
    inactive_cells = []
    offshore_boundary_cells = []
    onshore_boundary_cells = []
    #  surface_boundary_cells = []

    # add inactive cells
    for i in range(int(pars.ncol*pars.offshore_proportion)):
        for j in range(pars.nrow):
            for k in range(0, int((pars.Lz-pars.sea_level)/delv)):
                inactive_cells.append([k, j, i])
        
    # add cells on ends of domain
    for k in range(pars.nlay):
        for j in range(pars.nrow):
            if k >= np.floor((pars.Lz-pars.sea_level)/delv): # if the cell is below sea level
                offshore_boundary_cells.append([k, j, 0])

            if k >= np.floor((pars.Lz-pars.sea_level-pars.h_b)/delv):
                onshore_boundary_cells.append([k, j, pars.ncol-1])

    # add the seafloor
    for i in range(int(pars.ncol*pars.offshore_proportion)):
         for j in range(pars.nrow):
            offshore_boundary_cells.append([int((pars.Lz-pars.sea_level)/delv), j, i])

    # # add the recharge surface 
    # for i in range(int(offshore_proportion*ncol), ncol):
    #     for j in range(nrow):
    #         surface_boundary_cells.append([0, j, i])

    # create ibound array
    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)
    for cell in inactive_cells:
        ibound[cell[0], cell[1], cell[2]] = 0
    for cell in onshore_boundary_cells+offshore_boundary_cells:
        ibound[cell[0], cell[1], cell[2]] = -1

    # define starting heads
    strt = pars.sea_level*np.ones((pars.nlay, pars.nrow, pars.ncol))

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
        )

    # create solver package
    pcg = flopy.modflow.ModflowPcg(
            swt, 
            hclose=1.0e-5, 
            npcond=1, 
            mxiter=500
        )
    

    oc_spd = {} 
    for kstp in range(0, int(pars.perlen/pars.dt),pars.frequency):
            oc_spd[(0, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # Set up source sink data
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []

    # define onshore boundary data
    for cell in onshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], pars.sea_level+pars.h_b, pars.sea_level+pars.h_b])

    # define offshore boundary data
    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], pars.sea_level, pars.sea_level])
        
    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1

    # define constant head package
    chd = flopy.modflow.ModflowChd(
            model=swt, 
            stress_period_data=chd_data, 
            ipakcb=ipakcb
        )

    # define recharge values
    # rech = {}
    # rech_sp1 = []
    # # rech_spd = np.zeros((pars.nrow, pars.ncol))
    # # rech_spd[:, int(pars.ncol*pars.offshore_proportion):] = pars.W_net*np.ones((pars.nrow, int(pars.ncol-pars.ncol*pars.offshore_proportion)))

    # for j in range(pars.nrow):
    #     for i in range(int(pars.ncol*pars.offshore_proportion), ncol):


    # rech[0] = rech_sp1
    # create recharge package
    rch = flopy.modflow.ModflowRch(
            model=swt,
            rech=pars.W_net, #rech,
            ipakcb=ipakcb
        )

    # set starting concentrations
    sconc = 35.0*np.ones((pars.nlay, pars.nrow, pars.ncol))

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
    adv = flopy.mt3d.Mt3dAdv(
            model=swt, 
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
            percel=0.5
        )

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
                        pars.nrow*pars.ncol))

    # define source sink mixing package
    ssm = flopy.mt3d.Mt3dSsm(
            model=swt,
            stress_period_data=ssm_data, 
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
        os.mkdir(ws)

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