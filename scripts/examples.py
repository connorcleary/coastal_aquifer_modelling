from coastal_aquifer_model import DrainedCoastalAquifer

def steady_example():
    # example of setting up a steady state model
    # you just need to make sure you change the executable path :)
    model = DrainedCoastalAquifer("steady", exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    # set up model structure, click into function to see defaults
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
    model.plot_results()

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
    # set up sea level rise, again see in the function for defaults
    model.set_sea_level_rise(sea_level_rise=1)

    model.build_model()
    model.run_model()
    model.extract_results()
    model.plot_sea_level_rise_results()
    model.plot_drain_discharge_and_salinity()

def flux_controlled_example():

    head_controlled_model = DrainedCoastalAquifer("steady")
    inland_influx = head_controlled_model.get_final_inland_influx()

    steady_model = DrainedCoastalAquifer("steady_fixed_flux", exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    # set up steady_model structure
    steady_model.set_geometries(nlay=30, ncol=100)
    # set up temporal parameters
    steady_model.set_temporal_params()
    # set boundary conditions
    steady_model.set_boundary_conditions(h_w=0.5, q_b=inland_influx)
    # set hydraulic properties
    steady_model.set_hydraulic_properties()
    # set transport properties
    steady_model.set_transport_properties()
    # set up drain
    steady_model.set_up_drain(x_w=100, z_w=2)

    steady_model.build_model()
    steady_model.run_model()
    steady_model.extract_results()
    steady_model.plot_evolutions()
    steady_model.plot_results()

    slr_model = DrainedCoastalAquifer("slr_fixed_flux", exe_path=r"C:\Users\conno\Desktop\swt_v4_00_05\swt_v4_00_05\exe\swt_v4.exe")
    # set up slr_model structure
    slr_model.set_geometries(nlay=30, ncol=100)
    # set up temporal parameters
    slr_model.set_temporal_params()
    # set boundary conditions
    slr_model.set_boundary_conditions(h_w=0.5, q_b=inland_influx)
    # set hydraulic properties
    slr_model.set_hydraulic_properties()
    # set transport properties
    slr_model.set_transport_properties()
    # set up drain
    slr_model.set_up_drain(x_w=100, z_w=2)
    # set up sea level rise
    slr_model.set_sea_level_rise(sea_level_rise=1)

    slr_model.build_model()
    slr_model.run_model()
    slr_model.extract_results()
    slr_model.plot_sea_level_rise_results()
    slr_model.plot_drain_discharge_and_salinity()


def main():
    # steady_example()
    # sea_level_rise_example()
    flux_controlled_example()



if __name__ == "__main__":
    main()