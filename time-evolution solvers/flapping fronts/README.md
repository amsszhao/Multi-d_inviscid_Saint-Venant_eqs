## Description

In general, one should first run the python file XXX.py to generate raw data files and then run XXX.m file to generate movies. In all XXX.py files, the froude number is set to be 6 and the minimum fluid height is set to be 0.28 and the channel width and the y-boundary condition are varied.

Movies are downloadable from https://zenodo.org/record/8285097.

#### roll_width_point15.py -> roll_width_point15.m

With wall boundary condition, the width of the channel is set to be 0.15. No flapping front is seen at the end. Channel roll waves are stable.


#### roll_width_point16.py -> roll_width_point16.m

With wall boundary condition, the width of the channel is set to be 0.16. No flapping front is seen at the end. Channel roll waves are stable.


#### roll_width_point17.py -> roll_width_point17.m

With wall boundary condition, the width of the channel is set to be 0.17. No flapping front is seen at the end. Channel roll waves are stable.

#### roll_width_point18.py -> roll_width_point18.m

With wall boundary condition, the width of the channel is set to be 0.18. Wave fronts start to flap after a while. The flapping waves are persistent and do not become chaotic.

#### roll_width_point18_refined.py -> roll_width_point18_refined.m

With wall boundary condition, the width of the channel is set to be 0.18. Wave fronts start to flap after a while. The flapping waves are persistent and do not become chaotic. Finer meshgrid is used compared with roll_width_point18.py. Raw data files are used to generate figures in the paper.

#### roll_width_point18_periodic.py -> roll_width_point18_periodic.m

With periodic y-boundary condition, the width of the channel is set to be 0.18. No flapping front is seen at the end. Channel roll waves are stable.

#### roll_width_point36_periodic.py -> roll_width_point36_periodic.m

With periodic y-boundary condition, the width of the channel is set to be 0.36. Wave fronts start to flap after a while. The flapping waves are persistent and do not become chaotic.

#### roll_width_point2.py -> roll_width_point2.m

With wall boundary condition, the width of the channel is set to be 0.2. Wave fronts start to flap after a while. The flapping waves are also unstable, transitioning to chaotic flow at the end.

#### roll_width_1.py -> roll_width_1.m

With wall boundary condition, the width of the channel is set to be 1. Wave fronts start to flap after a while. The flapping waves are also unstable, transitioning to chaotic flow at the end.




