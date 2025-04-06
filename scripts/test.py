from coastal_aquifer_model import DrainedCoastalAquifer

def steady_example():
    model = DrainedCoastalAquifer("steady", exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    # set up model structure
    model.set_geometries(nlay=30, ncol=100)
    # set up temporal parameters
    model.set_temporal_params()
    # set boundary conditions
    model.set_boundary_conditions(h_w=0.5, h_b=1)
    # set hydraulic properties
    model.set_hydraulic_properties()
    # set transport properties
    model.set_transport_properties()
    # set up drain
    model.set_up_drain(x_w=100, z_w=2)

    model.build_model()
    model.run_model()
    model.extract_results()
    model.plot_evolutions()
    model.plot_boundary_concentration()
    model.plot_salinity()

def sea_level_rise_example():

    # initiate model
    model = DrainedCoastalAquifer("sea_level_rise", exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    # set up model structure
    model.set_geometries(nlay=30, ncol=100)
    # set up temporal parameters
    model.set_temporal_params()
    # set boundary conditions
    model.set_boundary_conditions(h_w=0.5, h_b=1)
    # set hydraulic properties
    model.set_hydraulic_properties()
    # set transport properties
    model.set_transport_properties()
    # set up drain
    model.set_up_drain(x_w=100, z_w=2)
    # set up sea level rise
    model.set_sea_level_rise(sea_level_rise=1)

    model.build_model()
    model.run_model()
    model.extract_results()
    model.plot_sea_level_rise_results()
    model.plot_drain_discharge_and_salinity()


def main():
    steady_example()
    sea_level_rise_example()



if __name__ == "__main__":
    main()