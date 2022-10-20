import os
import pickle

class ModelParameters:
    """
        Class to store model parameters

        Inputs:
            name: The name of the model/scenario/realization
            Lx: Length of the aquifer model [m]
            Ly: Alongshore dimension - set as 1 for the 2D case [m]
            Lz: Depth of the model [m]
            offshore_proportion: portion of the aquifer which is submarine 
            sea_level: datum for sea level (gives space for water table mound) [m]
            ncol: number of columns
            nrow: number of rows
            nlay: number of layers
            K: Horizontal hydraulic conductivity [m/day]
            anis: ratio between K and vertical hydraulic conductivity
            sy: specific yield 
            ss: specific storage [m^-1]
            n: porosity
            alpha_L: longitudinal dispersivity
            alpha_anisT: ratio between longitudinal and transverse dispersivity
            alpha_anisV: ratio between longitudinal and vertical dispersivity
                !!! I've set this as 0.1, but which werner says is the transverse ratio:
                    but there is no mention of vertical dispersivity !!!
            diff: molecular diffusion coefficient
            perlen: simulation period length [days]
            dt: timestep [days]
            h_b: inland boundary head
            W_net: recharge [m\day]
            rho_f: fresh water density [kg/m^3]
            rho_s: salt water density [kg/m^3]
            exe_path: path to the seawat executable
            frequency: frequency of timesteps save for results
    """
    
    def __init__(self, name="none", Lx=200, Ly=1, Lz=5.5, offshore_proportion=0.025, 
                sea_level=5, ncol=400, nrow=1, nlay=110, K=10, anis=1, sy=0.24, 
                ss=1e-5, n=0.3, alpha_L=1, alpha_anisT=0.1, alpha_anisV=0.01, 
                diff=8.64e-5, perlen=1e5, dt=1e2, h_b=0, W_net=0.00285,
                rho_f=1000, rho_s=1025, exe_file=r".\exe_path.txt", frequency=1, 
                x_w=0, Lx_w=10, Ly_w=1, z_w=0, h_w=0):


        self.name=name
        self.Lx=Lx
        self.Ly=Ly
        self.Lz=Lz
        self.offshore_proportion=offshore_proportion
        self.sea_level=sea_level
        self.ncol=ncol
        self.nrow=nrow
        self.nlay=nlay
        self.K=K
        self.anis=anis 
        self.sy=sy
        self.ss=ss
        self.n=n
        self.alpha_L=alpha_L
        self.alpha_anisT=alpha_anisT
        self.alpha_anisV=alpha_anisV
        self.diff=diff
        self.perlen=perlen
        self.dt=dt
        self.h_b=h_b
        self.W_net=W_net
        self.rho_f=rho_f
        self.rho_s=rho_s
        with open(exe_file, 'r') as f: self.exe_path=f.readline()
        self.frequency=frequency
        self.nstp=perlen/dt
        self.x_w=x_w
        self.Lx_w=Lx_w
        self.Ly_w=Ly_w
        self.z_w=z_w
        self.h_w=h_w

        self.save_parameters()

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d

    def save_parameters(self):
        """
            Save object
        """
        model_ws = f".\\model_files\\{self.name}"
        if not os.path.exists(model_ws):
            os.makedirs(model_ws)

        f = open(f"{model_ws}\\pars", 'wb')
        pickle.dump(self, f)
        f.close()

def load_parameters(name):

    model_ws = f".\\model_files\\{name}"
    f = open(f"{model_ws}\\pars", 'rb')
    temp = pickle.load(f)
    f.close()          
    return temp