## Description 

### point_sub.m

This is the main file. It realizes the computation of the Evans-Lopatinsky determinant of a discontinuous hydraulic shock with Froude number F and right fluid height HR.

### compute_V_n.m compute_V0_n.m compute_V.m compute_V0.m interior_ode_nonradial.m interior0_ode_nonradial.m interior_ode.m interior0_ode.m static_constant_v.m

These are supportive files of the main file point_sub.m.

### parameterspace.m

This generates a set of parameters of discontinuous hydraulic shocks for the investigation of spectral stability.

### extract_sub.m

This function helps to coordinate parallel computation and will extract useful stability info.

### job1.m ... job15.m

These are batch jobs sent to execute on supercomputers at Indiana University.

