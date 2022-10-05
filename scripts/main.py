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
    swt = cam.build_steady_model(pars)
    cam.run_model(swt)
    concentration, head, qx, qy, qz = cam.extract_results(name)
    results.plot_results(name)
    results.plot_evolutions(name)
    results.save_metrics(name, fraction=0.05)
    results.plot_boundary_concentration(name)


def main():
    create_run_plot_model("elongated", Lx=2000)

if __name__=="__main__":
    main()