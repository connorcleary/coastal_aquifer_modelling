from coastal_aquifer_model import DrainedCoastalAquifer


def create_run_plot_model(name):
    """
        Create, run and plot a new scenario

        Inputs:
            name: name for the model
            **kwargs: any key word arguments to be passed to the ModelParameters class
        Outputs:
            None
    """

    model = DrainedCoastalAquifer(name, exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    model.set_geometries()
    model.set_temporal_params()
    model.set_boundary_conditions()
    model.set_hydraulic_properties()
    model.set_transport_properties()
    model.set_up_drain(x_w=100, z_w=2)
    model.set_sea_level_rise(sea_level_rise=1)

    model.build_model()
    model.run_model()
    model.extract_results()
    model.plot_evolutions()
    model.save_metrics(fraction=0.005)
    model.plot_boundary_concentration()
    model.plot_sea_level_rise_results()
    model.plot_drain_discharge_and_salinity()


def main():
    create_run_plot_model('test_w_slr_drn')


if __name__ == "__main__":
    main()