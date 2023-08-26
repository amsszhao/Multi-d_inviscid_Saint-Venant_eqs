## Description

One should first run python files to generate raw date files and then run matlab files to generate movies.

Movies are downloadable from https://zenodo.org/record/8285097.

#### dam_break_F_equals_2_point_25.py -> dam_break_F_equals_2_point_25.m

Dam break initial data with left fluid height 1, right fluid height 0.7, and Froude number 2.25. The discontinuous hydraulic shock is 2d convectively unstable. Herringbone flows are seen at the end. 


#### flat_F_equals_2_point_25.py -> flat_F_equals_2_point_25.m

flat initial data with co-moving speed set to be that of discontinuous hydraulic shock with left limiting fluid height 1, right limiting fluid height 0.7, and  Froude number 2.25. Herringbone flows, parabola, and roll waves are seen at the end. 


#### dam_break_F_equals_2_point_13.py, dam_break_F_equals_2_point_14.py, dam_break_F_equals_2_point_15.py -> dam_break_comparison.m

Comparison of dam-break simulations with Froude number 2.13,2.14,2.15.

#### flat_F_equals_2_point_13.py, flat_F_equals_2_point_14.py, flat_F_equals_2_point_15.py -> flat_comparison.m

Comparison of flat simulations with Froude number 2.13,2.14,2.15.



