## Description 

### point_smooth.m

This is the main file. It realizes the computation of the Evans function of a smooth hydraulic shock with Froude number F and right fluid height HR.

### compute_V_n.m compute_V_nR.m compute_V.m compute_V0_n.m compute_V0_nR.m compute_V0.m interior_ode_nonradial.m interior0_ode_nonradial.m interior_ode.m interior0_ode.m static_constant_v.m

These are supportive files of the main file point_smooth.m.

### parameterspace.m

This generates a set of parameters of smooth hydraulic shocks for latter investigation of their spectral stability.

### extract_smooth.m

This function helps to coordinate parallel computation and will extract useful stability info.

### job1.m ... job23.m

These are batch jobs sent to execute on supercomputers at Indiana University.
