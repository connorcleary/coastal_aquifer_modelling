import coastal_aquifer_model as cam
from pars import ModelParameters, load_parameters
import results

def create_run_plot_model(name, **kwargs):
    """
        Create, run and plot a new scenario
        
        Inputs:
            name: name for the model
            **kwargs: any key word arguments to be passed to the ModelParameters class
        Outputs:
            None
    """

    pars = ModelParameters(name, **kwargs)
    swt = cam.build_model(pars)
    cam.run_model(swt)
    concentration, head, qx, qy, qz = cam.extract_results(name)
    results.plot_sea_level_rise_results(name)
    results.plot_evolutions(name)
    results.save_metrics(name, fraction=0.005)
    # results.plot_boundary_concentration(name)
    results.plot_drain_discharge_and_salinity(name)




def main():
    # create_run_plot_model("elongated", Lx=2000)
    # create_run_plot_model("case1_vka_sconc0_ncol4000_dt0.5_perlen1e3", h_b=0.3184, diff=0, alpha_L=0, W_net=0, ncol=4000, dt=0.5, perlen=1e3)
    create_run_plot_model("test_slr2", x_w=100, z_w=-2, h_w=0, h_b=0.5, nlay=40, drain_conductance=100000, Lz=6,
                          sea_level_rise=1, rise_type="linear", rise_length=365 * 100)

if __name__=="__main__":
    main()