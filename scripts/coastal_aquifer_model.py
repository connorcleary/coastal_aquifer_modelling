import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf
import pandas as pd
import matplotlib.pyplot as plt
import post_processing as proc

class DrainedCoastalAquifer():
    def __init__(self, name, exe_path):
        self.name = name
        self.drain = False
        self.exe_path=exe_path
        self._hydraulic_properties_set = False
        self._sea_level_rise_set = False

    def set_geometries(self, Lx=200, Ly=1, Lz=6.5, offshore_boundary_type='vertical', offshore_proportion=0.025,
                sea_level=5, ncol=400, nrow=1, nlay=110):
        """

        :param Lx: length of the aquifer [m]
        :param Ly: width of the aquifer [m]
        :param Lz: depth of the aquifer [m]
        :param offshore_boundary_type: "werner" or "vertical"
        :param offshore_proportion: if werner, fraction of the aquifer which is under the seafloor
        :param sea_level: sea level relative to the base of the aquifer [m]
        :param ncol:
        :param nrow:
        :param nlay:
        """
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        assert offshore_boundary_type in ['vertical', 'werner']
        self.offshore_boundary_type = offshore_boundary_type
        if self.offshore_boundary_type == 'vertical':
            self.offshore_proportion = 0
        else:
            self.offshore_proportion = offshore_proportion
        self.sea_level = sea_level
        self.ncol = ncol
        self.nrow = nrow
        self.nlay = nlay
        self.delr = Lx / ncol
        self.delc = Ly / nrow
        self.delv = Lz / nlay
    
    def set_hydraulic_properties(self, K=10, anis=1, sy=0.24,  ss=1e-5, n=0.3):
        """

        :param K: aquifer hydraulic conductivity [m/day]
        :param anis: vertical anisotropy of hydraulic conductivity
        :param sy: specific yield
        :param ss: specific storage
        :param n: porosity
        """
        self.K = K
        self.anis = anis
        self.sy = sy
        self.ss = ss
        self.n = n
    
    def set_transport_properties(self, alpha_L=1, alpha_anisT=0.1, alpha_anisV=0.01, 
                diff=8.64e-5, rho_f=1000, rho_s=1025):
        """

        :param alpha_L: longitudinal dispersivity
        :param alpha_anisT: ratio of transverse to longitudinal dispersivity
        :param alpha_anisV: ratio of vertical to longitudinal dispersivity
        :param diff: diffussion coefficient
        :param rho_f: density of fresh water
        :param rho_s: density of sea water
        """
        self.alpha_L = alpha_L
        self.alpha_anisT = alpha_anisT
        self.alpha_anisV = alpha_anisV
        self.diff = diff
        self.rho_f = rho_f
        self.rho_s = rho_s
        self._hydraulic_properties_set = True
        
    def set_temporal_params(self, perlen=1e5, dt=1e2, frequency=1):
        """

        :param perlen: period length for first (steady state period) [days]
        :param dt: time step [days]
        :param frequency: frequency to save results
        """
        self.perlen = perlen
        self.dt = dt
        self.frequency = frequency
        
    def set_up_drain(self, x_w=0, Lx_w=10, Ly_w=1, z_w=2, wetland_as_drain=True, drain_conductance=None):
        """

        :param x_w: distance from the offshore boundary to wetland [m]
        :param Lx_w: length of wetland [m]
        :param Ly_w: width of wetland [m]
        :param z_w: depth of wetland below land surface [m]
        :param wetland_as_drain: bool
        :param drain_conductance:
        """
        assert self._hydraulic_properties_set, "Please set hydraulic properties first"
        self.x_w = x_w
        self.Lx_w = Lx_w
        self.Ly_w = Ly_w
        self.z_w = z_w
        self.wetland_as_drain = wetland_as_drain
        if self.wetland_as_drain: 
            if drain_conductance == None:
                self.drain_conductance = self.K
            else:
                self.drain_conductance = drain_conductance
        self.drain = True
        
    def set_boundary_conditions(self, h_b=0, W_net=0.00285, h_w=0):
        """

        :param h_b: inland boundary head (above steady state sea-level)
        :param W_net: distributed recharge
        :param h_w: wetland head (above steady state sea-level
        """
        self.h_b = h_b
        self.W_net = W_net
        self.h_w = h_w
        
    def set_sea_level_rise(self, sea_level_rise=1, rise_type='linear', rise_length=100*365, rise_time_series=None, initial_conditions=None):
        """

        :param sea_level_rise: sea level rise amout [m]
        :param rise_type: type of rise ("linear" or "step")
        :param rise_length: length of sea-level rise period [days]
        :param rise_time_series:
        :param initial_conditions:
        """
        assert rise_type in ["linear", "step"]
        self.sea_level_rise = sea_level_rise
        self.rise_type = rise_type
        self.rise_length = rise_length
        self.rise_time_series = rise_time_series
        self.initial_conditions = initial_conditions
        self._sea_level_rise_set = True
        
    def _create_temporal_discretization(self):
        """
        Funciton to make the tdis period data
        :param self:
        :param sea_level_rise: sea level rise [m]
        :param rise_length: number of days for sea level rise period [days]
        :param rise_type: type of rise "linear" or "step"
        :param rise_rate:
        :param initial_conditions:
        :return:
        """
        perlen = []
        nstep = []
        steady = []
        if self.initial_conditions is None:
            perlen.append(self.perlen)
            nstep.append(self.perlen/self.dt)
            steady.append(True)
            if self.sea_level_rise == 0:
                return perlen, nstep, steady, len(perlen)
        else:
            raise NotImplementedError
    
        if self.rise_type == 'step':
            perlen.append(self.rise_length)
            nstep.append(np.round(self.rise_length/self.dt))
            steady.append(False)
            return perlen, nstep, steady, len(perlen)
        if self.rise_type == 'linear':
            for i in range(int(self.rise_length/self.dt)):
                perlen.append(self.dt)
                nstep.append(1)
                steady.append(False)
            return perlen, nstep, steady, len(perlen)
        
    def _create_cell_groups(self):
        """
        Create lists of cell groups
        :param self:
        :param delr:
        :param delv:
        :param delc:
        :return:
        """
        # define cell groups
        inactive_cells = []
        offshore_boundary_cells = []
        onshore_boundary_cells = []
        #  surface_boundary_cells = []
        wetland_cells = []
    
        # add inactive cells
        if self.offshore_boundary_type == "werner":
            for i in range(int(self.ncol * self.offshore_proportion)):
                for j in range(self.nrow):
                    for k in range(0, int((self.Lz - self.sea_level) / self.delv)):
                        inactive_cells.append([k, j, i])
    
        # add cells on ends of domain
        for k in range(self.nlay):
            for j in range(self.nrow):
                if k >= np.floor((self.Lz - self.sea_level) / self.delv):  # if the cell is below sea level
                    offshore_boundary_cells.append([k, j, 0])
    
                if k >= np.floor((self.Lz - self.sea_level - self.h_b) / self.delv):
                    onshore_boundary_cells.append([k, j, self.ncol - 1])
    
        # add the seafloor
        if self.offshore_boundary_type == "werner":
            for i in range(int(self.ncol * self.offshore_proportion)):
                for j in range(self.nrow):
                    offshore_boundary_cells.append([int((self.Lz - self.sea_level) / self.delv), j, i])
    
        # add wetland cells
        if self.x_w != 0:
            for j in range(self.nrow):
                for i in range(int((self.x_w) / self.delr),
                               int((self.x_w + self.Lx_w) / self.delr)):
                    for k in range(int((self.Lz - self.sea_level - self.h_w) / self.delv) + 1,
                                   int((self.Lz - self.sea_level - self.z_w) / self.delv)):
                        wetland_cells.append([k, j, i])
                    for k in range(0, int((self.Lz - self.sea_level - self.h_w) / self.delv) + 1):
                        inactive_cells.append([k, j, i])
    
        # create ibound array
        ibound = np.ones((self.nlay, self.nrow, self.ncol), dtype=np.int32)
        for cell in inactive_cells:
            ibound[cell[0], cell[1], cell[2]] = 0
        else:
            for cell in onshore_boundary_cells + offshore_boundary_cells:
                ibound[cell[0], cell[1], cell[2]] = -1
                
        return inactive_cells, offshore_boundary_cells, onshore_boundary_cells, wetland_cells, ibound
    
    def _create_stress_period_data(self, wetland_cells, onshore_boundary_cells, offshore_boundary_cells):
        """
        Create stress period data for all packages
    
        :param self:
        :param initial_conditions:
        :param wetland_cells:
        :param onshore_boundary_cells:
        :param offshore_boundary_cells:
        :return:
        """
        drn_spd = {per: [] for per in range(self.nper)}
        chd_spd = {per: [] for per in range(self.nper)}
        ssm_spd = {per: [] for per in range(self.nper)}
        oc_spd = {}
        oc_spd = {}
        itype = flopy.mt3d.Mt3dSsm.itype_dict()
        for kstp in range(0, int(self.perlen / self.dt), self.frequency):
            oc_spd[(0, kstp)] = ["save head", "save budget"]
    
        for cell in wetland_cells:
            for per in range(self.nper):
                drn_spd[per].append([cell[0], cell[1], cell[2], self.sea_level+self.h_w, self.drain_conductance])
                ssm_spd[per].append([cell[0], cell[1], cell[2], 0.0, itype["DRN"]])
    
        for cell in onshore_boundary_cells:
            for per in range(self.nper):
                ssm_spd[per].append([cell[0], cell[1], cell[2], 0, itype["BAS6"]])
                chd_spd[per].append([cell[0], cell[1], cell[2], self.sea_level+self.h_b, self.sea_level+self.h_b])
    
        if self.initial_conditions is None:
            for cell in offshore_boundary_cells:
                ssm_spd[0].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
                chd_spd[0].append([cell[0], cell[1], cell[2], self.sea_level, self.sea_level])
                for kstp in range(0, int(self.perlen / self.dt), self.frequency):
                    oc_spd[(0, kstp)] = ["save head", "save budget"]
    
        if self.rise_type == 'step':
            for cell in offshore_boundary_cells: #+ _get_extra_chd_cells(self, sea_level_rise):
                ssm_spd[1].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
                chd_spd[1].append([cell[0], cell[1], cell[2], self.sea_level+self.sea_level_rise, self.sea_level+self.sea_level_rise])
                for kstp in range(0, int(self.rise_length / self.dt), self.frequency):
                    oc_spd[(1, kstp)] = ["save head", "save budget"]
    
        if self.rise_type == 'linear':
            for per in range(1, self.nper):
                sea_level_rise_increment = per/self.nper*self.sea_level_rise
                for cell in offshore_boundary_cells: #+ _get_extra_chd_cells(self, sea_level_rise_increment):
                    ssm_spd[per].append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
                    chd_spd[per].append([cell[0], cell[1], cell[2], self.sea_level + sea_level_rise_increment,
                                       self.sea_level + sea_level_rise_increment])
                if per%self.frequency == 0 or per == self.nper-1:
                    oc_spd[(per, 0)] = ["save head", "save budget"]
    
        return drn_spd, chd_spd, oc_spd, ssm_spd
    
    def build_model(self):
        '''
            A function to build a coastal aquifer model.
            :param self: parameters defining the geometry of the model
            :return:
        '''
        # create model workspace
        if not self._sea_level_rise_set:
            self.set_sea_level_rise(sea_level_rise=0)

        model_ws = f".\\model_files\\{self.name}"
        if not os.path.exists(model_ws):
            os.makedirs(model_ws)
    
        # create base seawat model
        swt = flopy.seawat.Seawat(self.name, model_ws=model_ws, exe_name=self.exe_path)
    
        # calculate cell dimension
        delr = self.Lx/self.ncol
        delc = self.Ly/self.nrow
        delv = self.Lz/self.nlay
    
        # define top and bottoms
        assert self.Lz-self.sea_level >= self.sea_level_rise, "Model top must be higher than sea-level"
        top = self.Lz
        botm = np.linspace(top-delv, 0, self.nlay)
    
        # something I've copied
        ipakcb = 53
    
        perlen, nstep, steady, nper = self._create_temporal_discretization()
        self.nper = nper
    
        # define discretization package
        dis = flopy.modflow.ModflowDis(
                model=swt,
                nlay=self.nlay,
                nrow=self.nrow,
                ncol=self.ncol,
                nper=nper,
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
    
        inactive_cells, offshore_boundary_cells, onshore_boundary_cells, wetland_cells, ibound = self._create_cell_groups()
    
        # define starting heads
        if self.initial_conditions is None:
            strt = self.Lz*np.ones((self.nlay, self.nrow, self.ncol))
        else:
            raise NotImplementedError
    
        # create basic package
        bas = flopy.modflow.ModflowBas(
                model=swt, 
                ibound=ibound, 
                strt=strt
            )
    
        # define layer types and wetting, this one I'm not sure about
        laytyp=np.ones(self.nlay)
        laytyp[0] = 1
        laywet=np.ones(self.nlay)
        laywet[0] = 1
    
        # create layer property flow package
        lpf = flopy.modflow.ModflowLpf(
                swt, 
                hk=self.K, 
                vka=self.anis, 
                ipakcb=ipakcb, 
                laytyp=laytyp, 
                laywet=laywet,
                ss=self.ss, # not sure about these ones
                sy=self.sy,
                layvka=1,
            )
    
        # create solver package
        pcg = flopy.modflow.ModflowPcg(
                swt, 
                hclose=1.0e-5, 
                npcond=1, 
                mxiter=500
            )
    
        drn_spd, chd_spd, oc_spd, ssm_spd = self._create_stress_period_data(wetland_cells, onshore_boundary_cells, offshore_boundary_cells)
    
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
                rech=self.W_net, #rech,
                ipakcb=ipakcb
            )
    
        # set starting concentrations
        if self.initial_conditions is None:
            sconc = 0.0*np.ones((self.nlay, self.nrow, self.ncol))
            sconc[:, :, 0] = 35.0
        else:
            raise NotImplementedError
    
        # define basic transport package
        btn = flopy.mt3d.Mt3dBtn(
                swt,
                nprs=-self.frequency,
                prsity=self.n,
                sconc=sconc,
                chkmas=False,
                nprobs=10,
                nprmas=10,
                dt0=self.dt
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
                al=self.alpha_L, 
                trpt=self.alpha_anisT, 
                trpv=self.alpha_anisV, 
                dmcoef=self.diff,
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
        if self.offshore_boundary_type == "werner":
            mxss = int(np.ceil(2*self.nlay*self.nrow +
                                self.nrow*self.ncol*self.offshore_proportion+1 +
                                self.nrow*self.ncol+self.Lx_w*self.nrow*(self.x_w!=0)))
        else:
            mxss = int(np.ceil(2 * self.nlay * self.nrow +
                               self.nrow * self.ncol + self.Lx_w * self.nrow * (self.x_w != 0)))

    
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
                firstdt=self.dt,
            )
    
        # write input
        swt.write_input() 
    
        self.swt = swt
    
    def run_model(self):
        """
            A function to run the seawat model
    
            Inputs: 
                swt: model object
            Outputs:
                None
        """
        self.swt.write_input()
        success, buff = self.swt.run_model(silent=False, report=True)
        if not success:
            raise Exception("SEAWAT did not terminate normally.")
    
    
    def extract_results(self):
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
        name = self.name
        model_ws = f".\\model_files\\{name}"
        nstp = self.perlen/self.dt
    
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
            if self.nrow > 1:
                qy[t] = cbbobj.get_data(text="flow front face", totim=times[t])[0]
            qz[t] = cbbobj.get_data(text="flow lower face", totim=times[t])[0]
    
        self.save_results(concentration, head, qx, qy, qz)
    
    
    def save_results(self, concentration, head, qx, qy, qz):
        """
            Save extracted results to a .npy file
    
            Inputs:
                name: model name
                concentration, head etc. : numpy arrays of model outputs
            Outputs:
                None
        """
        name = self.name
        ws = os.path.join(f'.\\results\\{name}')
        if not os.path.exists(ws):
            os.makedirs(ws)
    
        with open(os.path.join(ws, f"qx.npy"), 'wb') as f: np.save(f, np.array(qx))
        with open(os.path.join(ws, f"qy.npy"), 'wb') as f: np.save(f, np.array(qy))
        with open(os.path.join(ws, f"qz.npy"), 'wb') as f: np.save(f, np.array(qz))
        with open(os.path.join(ws, f"head.npy"), 'wb') as f: np.save(f, np.array(head))
        with open(os.path.join(ws, f"concentration.npy"), 'wb') as f: np.save(f, np.array(concentration))
    
    
    def load_results(self):
        """
            Load extracted results from .npy files
    
            Inputs:
                name: name of the model
            Outputs:
                concentration, head... : numpy matrices of results
        """
        name = self.name
        ws = os.path.join(f'.\\results\\{name}')
        if not os.path.exists(os.path.join(ws, f"qx.npy")):
            self.extract_results()
    
        with open(os.path.join(ws, f"qx.npy"), 'rb') as f: qx = np.load(f, allow_pickle=True)
        with open(os.path.join(ws, f"qy.npy"), 'rb') as f: qy = np.load(f, allow_pickle=True)
        with open(os.path.join(ws, f"qz.npy"), 'rb') as f: qz = np.load(f, allow_pickle=True)
        with open(os.path.join(ws, f"head.npy"), 'rb') as f: head = np.load(f, allow_pickle=True)
        with open(os.path.join(ws, f"concentration.npy"), 'rb') as f: concentration = np.load(f, allow_pickle=True)
    
        return concentration, head, qx, qy, qz,
    
    def get_drain_results(self):
        ws = os.path.join(f'.\\results\\{self.name}')
        if not os.path.exists(os.path.join(ws, f"drain.npz")):
            self._extract_drain_results()
    
        results = np.load(os.path.join(ws, f"drain.npz"))
    
        return results
    
    def _extract_drain_results(self):
    
        name = self.name
        model_ws = f".\\model_files\\{name}"
    
        # open binary files
        cbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
        ucnobject = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    
        times = cbobj.get_times()
        drain_discharge = []
        drain_conc = []
        _, _, _, wetland_cells, _ = self._create_cell_groups()
    
        for t in times:
    
            temp_discharge = cbobj.get_data(text="DRAINS", totim=t)[0]
            drain_discharge.append(np.sum(temp_discharge['q']))
    
            conc = ucnobject.get_data(totim=t)
            conc_times_discharge = 0
            for flow in temp_discharge:
                shape = (self.nlay, self.nrow, self.ncol)
                cell = np.unravel_index(flow['node'], shape)
                conc_times_discharge += flow['q']*conc[cell[0], cell[1], cell[2]]
    
            drain_conc.append(conc_times_discharge/drain_discharge[-1])
    
            np.savez(f".\\results\\{name}\\drain.npz", discharge=drain_discharge, conc=drain_conc, times=times)

    def _find_mixing_volume(self, conc, fraction=0.05):
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
        for i in range(self.nlay):
            for j in range(self.ncol):
                if 35. * fraction <= conc[i, j] <= 35 * (1 - fraction):
                    mix += 1
        mixing_zone = self.Ly * (self.Lx / self.ncol) * (self.Lz / self.nlay) * mix

        return mixing_zone

    def _find_toe_position(self, conc, fraction=0.05):
        """
            Find the toe position

            Inputs:
                conc: 2D concentation array
                fraction: what fraction of salt water to consider
                pars: ModelParameters object
            Outputs:
                toe: position of the toe
        """

        for i in range(self.ncol - 1, -1, -1):
            if conc[-1][i] > 35 * fraction:
                return (i - self.ncol * self.offshore_proportion) * (self.Lx / self.ncol)

    def _find_mixing_centroid(self, conc, fraction=0.05):
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

        for i in range(self.nlay):
            for j in range(self.ncol):
                if 35. * fraction <= conc[i, j] <= 35 * (1 - fraction):
                    mix_n += 1
                    x_tot += (j - self.ncol * self.offshore_proportion) * (self.Lx / self.ncol)

        if mix_n != 0:
            centroid = x_tot / mix_n
        else:
            centroid = 0

        return centroid

    def _find_boundary_fluxes(self, conc, qx, qz, fraction=0.05):
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

        delv = self.delv

        offshore_inflow_s = 0
        offshore_inflow_f = 0
        offshore_outflow_s = 0
        offshore_outflow_f = 0
        onshore_inflow_s = 0
        onshore_inflow_f = 0
        onshore_outflow_s = 0
        onshore_outflow_f = 0

        # look at cells on the ends of the domain
        for k in range(self.nlay):
            # offshore cells
            if k >= np.floor((self.Lz - self.sea_level) / delv):  # if the cell is below sea level
                # if saline
                if conc[k, 0] >= 35 * fraction:
                    # if influx
                    if qx[k, 0] > 0:
                        offshore_inflow_s += qx[k, 0]
                    # if outflux
                    else:
                        offshore_outflow_s -= qx[k, 0]
                # if fresh
                else:
                    # if influx
                    if qx[k, 0] > 0:
                        offshore_inflow_f += qx[k, 0]
                    # if outflux
                    else:
                        offshore_outflow_f -= qx[k, 0]

            # onshore cells
            if k >= np.floor((self.Lz - self.sea_level - self.h_b) / delv):
                # if saline
                if conc[k, -2] >= 35 * fraction:
                    # if influx
                    if qx[k, -2] < 0:
                        onshore_inflow_s -= qx[k, -2]
                    # if outflux
                    else:
                        onshore_outflow_s += qx[k, -2]
                # if fresh
                else:
                    # if influx
                    if qx[k, -2] < 0:
                        onshore_inflow_f -= qx[k, -2]
                    # if outflux
                    else:
                        onshore_outflow_f += qx[k, -2]

        # look at seafloor cells
        for i in range(int(self.ncol * self.offshore_proportion)):
            # if saline
            if conc[int((self.Lz - self.sea_level) / delv), i] >= 35 * fraction:
                # if influx
                if qz[int((self.Lz - self.sea_level) / delv), i] > 0:
                    offshore_inflow_s += qz[int((self.Lz - self.sea_level) / delv), i]
                # if outflux
                else:
                    offshore_outflow_s -= qz[int((self.Lz - self.sea_level) / delv), i]
            # if fresh
            else:
                # if influx
                if qz[int((self.Lz - self.sea_level) / delv), i] > 0:
                    offshore_inflow_f += qz[int((self.Lz - self.sea_level) / delv), i]
                # if outflux
                else:
                    offshore_outflow_f -= qz[int((self.Lz - self.sea_level) / delv), i]

        return offshore_inflow_s, offshore_inflow_f, offshore_outflow_s, offshore_outflow_f, \
            onshore_inflow_s, onshore_inflow_f, onshore_outflow_s, onshore_outflow_f

    def _find_wetland_flux(self, conc, qx, qz):

        conc_flux = 0

        delv = self.Lz / self.nlay
        delr = self.Lx / self.ncol
        wetland_base = [[int((self.Lz - self.sea_level - self.z_w) / delv), col] for col in
                        range(int((self.offshore_proportion * self.Lx + self.x_w) / delr),
                              int((self.offshore_proportion * self.Lx + self.x_w + self.Lx_w) / delr))]

        wetland_left = [[int(lay), int((self.offshore_proportion * self.Lx + self.x_w) / delr) - 1] for lay in
                        range(int((self.Lz - self.sea_level - self.h_w) / delv),
                              int((self.Lz - self.sea_level - self.z_w) / delv))]

        wetland_right = [[int(lay), int((self.offshore_proportion * self.Lx + self.x_w + self.Lx_w) / delr)] for lay in
                         range(int((self.Lz - self.sea_level - self.h_w) / delv),
                               int((self.Lz - self.sea_level - self.z_w) / delv))]

        for cell in wetland_base:
            conc_flux += np.max([0, -1 * qz[cell[0], cell[1]] * conc[cell[0], cell[1]]])

        for cell in wetland_left:
            conc_flux += np.max([0, qx[cell[0], cell[1]] * conc[cell[0], cell[1]]])

        for cell in wetland_right:
            conc_flux += np.max([0, -1 * qx[cell[0], cell[1]] * conc[cell[0], cell[1]]])

        return conc_flux

    def _find_mound(self, qx):
        """
            Find the position of the mound, based on the smallest horizontal flux
                on the water table

            Inputs:
                qx: flux array
                pars: ModelParameters object
        """
        min_flux_value = np.inf
        min_flux_position = 0

        for k in range(self.nlay):
            if float(np.max(qx[k, 1:])) != float(0):
                for i in range(self.ncol):
                    if float(qx[k, i]) != float(0) and np.abs(qx[k, i]) < min_flux_value:
                        min_flux_value = np.abs(qx[k, i])
                        min_flux_position = i
                break

        if self.W_net == 0:
            if self.h_b > 0:
                min_flux_position = self.ncol - 1
            else:
                min_flux_position = 0

        mound = (min_flux_position - self.ncol * self.offshore_proportion) * self.Lx / self.ncol

        return mound

    def plot_salinity(self, ax, conc, qx, qz, cmap="viridis", row=0, x_step=10, z_step=5, width=0.002, arrow_c="white"):

        x = np.linspace(-self.Lx * self.offshore_proportion, self.Lx - self.Lx * self.offshore_proportion, self.ncol)
        y = np.linspace(-self.sea_level, self.Lz - self.sea_level, self.nlay)

        concentration_array = conc[:, row, :]
        for i in range(self.nlay):
            for j in range(self.ncol):
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

    def plot_results(self, timestep=-1, row=0, return_axs=False, figsize=(12, 6),
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
        concentration, head, qx, qy, qz = self.load_results()

        f, axs = plt.subplots(2, 1, figsize=figsize)

        # set up x and y arrays, to be distance above sea level and distance onshore
        x = np.linspace(-self.Lx * self.offshore_proportion, self.Lx - self.Lx * self.offshore_proportion, self.ncol)
        y = np.linspace(-self.sea_level, self.Lz - self.sea_level, self.nlay)

        # select relevent slice in time and the alongshore direction, and set values above the water table as nan
        concentration_array = concentration[timestep, :, row, :]
        head_array = head[timestep, :, row, :] - self.sea_level * np.ones_like(head[timestep, :, row, :])
        for i in range(self.nlay):
            for j in range(self.ncol):
                if concentration_array[i, j] == np.float32(1.e30):
                    head_array[i, j] = np.nan
                    concentration_array[i, j] = np.nan

        # plot head colormesh
        headcm = axs[0].pcolormesh(x, y, np.flipud(head_array),
                                   cmap=cmap, vmax=np.nanmax(head_array[:, 1:]), vmin=np.min([0, self.h_b]))

        # plot head contours
        hc = axs[0].contour(x, y, np.flipud(head_array), colors=arrow_c,
                            levels=np.linspace(np.min([0, self.h_b]), np.nanmax(head_array[:, 1:]), 15))

        # label contours
        axs[0].clabel(hc, hc.levels, inline=True, fontsize=10, fmt=fmt)

        # plot concentration colormesh
        conccm = axs[1].pcolormesh(x, y, np.flipud(concentration_array),
                                   cmap=cmap, vmax=35, vmin=0)

        # plot arrows
        X, Y = np.meshgrid(x[::x_step], y[::z_step])
        axs[1].quiver(X, Y,
                      np.flipud(qx[timestep, ::z_step, row, ::x_step]),
                      np.flipud(qz[timestep, ::z_step, row, ::x_step]),
                      color=arrow_c, width=width)

        axs[0].set_aspect(aspect)
        axs[1].set_aspect(aspect)
        axs[0].set_title("Head")
        axs[1].set_title("Salinity")
        axs[0].set_ylabel("Height above sealevel (m)")
        axs[1].set_ylabel("Height above sealevel (m)")
        axs[1].set_xlabel("Distance onshore (m)")

        f.suptitle(f"Head and salinity distributions for {self.name}")

        headcb = plt.colorbar(headcm, shrink=1, ax=axs[0])
        conccb = plt.colorbar(conccm, shrink=1, ax=axs[1])
        headcb.ax.set_title('Head (m)', fontsize='small')
        conccb.ax.set_title('Salinity (kg/m^3)', fontsize='small')

        # return axs objects if necessary
        if return_axs:
            return axs
        else:
            ws = os.path.join(f'.\\figures\\{self.name}')
            if not os.path.exists(ws):
                os.makedirs(ws)
            plt.savefig(f"{ws}\\head_and_concentration", dpi=300)

    def plot_evolutions(self, row=0, fraction=0.01, return_axs=False, figsize=(18, 6), interval=20):
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
        concentration, head, qx, qy, qz = self.load_results()
        nstp = int(self.perlen / self.dt)
        # set as 1% saltwater

        # create arrays
        toe = np.zeros(int(nstp / interval))
        mixing_volume = np.zeros(int(nstp / interval))
        centroid = np.zeros(int(nstp / interval))

        # select every 20th step
        times = np.linspace(0, nstp - 1, int(nstp / interval), dtype=int)

        for tdx, t in enumerate(times):
            toe[tdx] = self._find_toe_position(concentration[t, :, row, :], fraction)
            mixing_volume[tdx] = self._find_mixing_volume(concentration[t, :, row, :], fraction)
            centroid[tdx] = self._find_mixing_centroid(concentration[t, :, row, :], fraction)

        f, axs = plt.subplots(3, 1, sharex=True, figsize=figsize)

        # plot lines
        axs[0].plot(self.dt * times, toe)
        axs[1].plot(self.dt * times, mixing_volume)
        axs[2].plot(self.dt * times, centroid)

        axs[0].set_ylabel("distance onshore (m)")
        axs[1].set_ylabel("volume (m^3)")
        axs[2].set_ylabel("distance onshore (m)")

        axs[2].set_xlabel("time (days)")

        axs[0].set_title("Toe position")
        axs[1].set_title("Mixing zone volume")
        axs[2].set_title("Mixing zone centroid")

        f.suptitle(f"Evolution of metrics for {self.name}")

        ws = os.path.join(f'.\\figures\\{self.name}')
        if not os.path.exists(ws):
            os.makedirs(ws)

        plt.savefig(f"{ws}\\metric_evolutions", dpi=300)

        if return_axs: return axs

    def plot_boundary_concentration(self, return_axs=False, figsize=(6, 6), row=0):
        """
            Plot the concentrations in the cells along the right hand boundary

            Inputs:
                name: model name
                return_axs: flag whether to return axes
                figsize: figure size in inches
                row: model row
        """

        concentration, head, qx, qy, qz = self.load_results()

        # find mound position
        mound = self._find_mound(qx[-1, :, row, :])
        mound_col = int(np.round((mound / self.Lx + self.offshore_proportion) * self.ncol))

        # truncate to just include layers below the inland boundary
        concentration_mound = np.copy(
            concentration[-1, int(self.nlay / self.Lz * (self.Lz - self.sea_level - self.h_b)):, row, mound_col])
        concentration_edge = np.copy(
            concentration[-1, int(self.nlay / self.Lz * (self.Lz - self.sea_level - self.h_b)):, row, -1])

        y = np.linspace(-self.sea_level, self.h_b, len(concentration_edge))

        inland_lay = int(self.nlay / self.Lz * (self.Lz - self.sea_level - self.h_b)) + 1
        # find isoclors at the mound column
        concentration_mound[-1] = 35
        lay_isoclor1 = np.atleast_1d(np.argmax(concentration_mound > 0.35))[0] + inland_lay
        lay_isoclor5 = np.atleast_1d(np.argmax(concentration_mound > 1.75))[0] + inland_lay
        lay_isoclor10 = np.atleast_1d(np.argmax(concentration_mound > 3.5))[0] + inland_lay

        f, ax = plt.subplots(figsize=figsize)
        ax.plot(np.flipud(concentration_edge), y)
        ax.set_ylabel("Distance above sea level (m)")
        ax.set_xlabel("Salinity (kg/m^3)")
        ax.set_title(f"Salinity at inland boundary for {self.name}")

        # plot isoclors
        delv = self.Lz / self.nlay
        ax.axhline(self.Lz - self.sea_level - (lay_isoclor1) * delv, c='b', alpha=0.5, zorder=-1, linestyle=':',
                   label=r"mound 1% isoclor")
        ax.axhline(self.Lz - self.sea_level - (lay_isoclor5) * delv, c='g', alpha=0.5, zorder=-1, linestyle='-',
                   label=r"mound 5% isoclor")
        ax.axhline(self.Lz - self.sea_level - (lay_isoclor10) * delv, c='r', alpha=0.5, zorder=-1, linestyle='--',
                   label=r"mound 10% isoclor")

        ax.set_xlim([-5, 40])
        ax.set_ylim([-self.sea_level - 0.5, self.Lz - self.sea_level + 0.5])
        ax.legend()

        ws = os.path.join(f'.\\figures\\{self.name}')
        if not os.path.exists(ws):
            os.makedirs(ws)

        plt.savefig(f"{ws}\inland_boundary_salinity", dpi=300)

        if return_axs: return ax

    def plot_sea_level_rise_results(self):
        """
        Plot series of snapshots of salinity cross sections
        :param name: 
        """
        assert self.sea_level_rise > 0, "Model needs to include sea level rise"
        f, axs = plt.subplots(4, 3, figsize=(10, 6))
        concentration, head, qx, qy, qz = self.load_results()
        axs = axs.flatten()
        if self.initial_conditions is None:
            i_init = int(self.perlen / self.dt / self.frequency)
            self.plot_salinity(axs[0], concentration[i_init], qx[i_init], qz[i_init])
            axs[0].set_title("Steady state")
        else:
            raise NotImplementedError("still need to set up using different initial conditions")

        for i in range(1, 11):
            step = int(self.rise_length / 10 / (self.dt * self.frequency))
            self.plot_salinity(axs[i], concentration[i_init + i * step], qx[i_init + i * step], qz[i_init + i * step])
            axs[i].set_title(f"{i * self.rise_length / 10 / 365:1.0f} years")

        axs[-1].set_axis_off()

        plt.tight_layout()
        ws = os.path.join(f'.\\figures\\{self.name}')
        if not os.path.exists(ws):
            os.makedirs(ws)

        plt.savefig(f"{ws}\concentration_snapshots.png", dpi=600)

    def plot_drain_discharge_and_salinity(self):
        """
        Plot discharge and salinity at the drain/wetland
        :param name: 
        """
        assert self.sea_level_rise > 0, "model needs to include sea level risee"
        f, axs = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
        if self.initial_conditions is not None:
            raise NotImplementedError

        results = self.get_drain_results()
        idx = results['times'] > results['times'][-1] - self.rise_length
        times = (results['times'][idx] - (results['times'][-1] - self.rise_length)) / 365
        discharges = results['discharge'][idx]
        concentrations = results['conc'][idx]

        axs[0].plot(times, -discharges)
        axs[0].set_ylabel("Discharge [m^2/day]")
        axs[1].plot(times, concentrations)
        axs[1].set_ylabel("Concentration [PSU]")
        axs[1].set_xlabel("Time [years]")

        ws = os.path.join(f'.\\figures\\{self.name}')
        if not os.path.exists(ws):
            os.makedirs(ws)

        plt.savefig(f"{ws}\drain_flow_and_concentration.png", dpi=600)

    def save_metrics(self, row=0, fraction=0.05):
        """
            Find and save metrics for the final timestep. Saves to text file

            Inputs:
                name: name of model
                row: row of model
                fraction: fraction to consider saline
            Outputs: 
                none
        """
        concentration, head, qx, qy, qz = self.load_results()

        offshore_inflow_s, offshore_inflow_f, offshore_outflow_s, offshore_outflow_f, \
            onshore_inflow_s, onshore_inflow_f, onshore_outflow_s, onshore_outflow_f = \
            self._find_boundary_fluxes(concentration[-1, :, row, :], \
                                      qx[-1, :, row, :], qz[-1, :, row, :], fraction=fraction)

        toe = self._find_toe_position(concentration[-1, :, row, :], fraction)
        mixing_volume = self._find_mixing_volume(concentration[-1, :, row, :], fraction)
        centroid = self._find_mixing_centroid(concentration[-1, :, row, :], fraction)

        mound = self._find_mound(qx[-1, :, row, :])

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

        units = ["m", "m", "m^3", "m", "m^3/s", "m^3/s", "m^3/s",
                 "m^3/s", "m^3/s", "m^3/s", "m^3/s", "m^3/s"]

        ws = os.path.join(f'.\\results\\{self.name}')
        if not os.path.exists(ws):
            os.makedirs(ws)

        with open(f"{ws}\\metrics_{fraction}.txt", "w") as f:

            f.write(f"Metrics for {self.name} with fraction={fraction} \n\n")
            for (metric, string, unit) in zip(metrics, strings, units):
                f.write(f"{string}: {metric} {unit}\n")

            if self.x_w != 0:
                f.write(
                    f"Concentration flux to wetland: {self._find_wetland_flux(concentration[-1, :, row, :], qx[-1, :, row, :], qz[-1, :, row, :])} kg m^-1 day^-1")

        concentration_b = concentration[-1, :, row, -1]
        # cell centres
        depths = np.linspace(self.Lz - self.sea_level - (self.Lz / (2 * self.nlay)),
                             -self.sea_level + (self.Lz / (2 * self.nlay)), self.nlay)
        lays = np.linspace(0, 109, 110)

        headers = ["layer", "masl", "C"]
        results = np.vstack((headers, np.stack((lays, depths, concentration_b)).T))
        np.savetxt(f"{ws}\\boundary_concentrations.csv", results, delimiter=",", fmt='% s')
