#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: dam break
==================================

Solve the 2D inclined shallow water equations:

.. math::
    h_t + (hu)_x + (hv)_y = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y = h-(u^2+v^2)^(1/2)*u \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y = -(u^2+v^2)^(1/2)*v .
    
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn

def qinit(state,HL=1.,HR=0.7,dam_break_location=16):
    c=(HR**(3/2) - HL)/(HR - HL)
    UL=HL**(1/2)-c
    UR=HR**(1/2)-c
    X, Y = state.p_centers
    perturbation=np.zeros(X.shape)
    for i in range(2000):
        for j in range(500):
            r2 =(X[i,j]-4.0)**2 + (Y[i,j]-2.5)**2
            if r2<=0.5:
                perturbation[i,j]=np.exp(-0.2/(0.5-r2))
    state.q[depth     ,:,:] = HL*(X<=dam_break_location) + HR*(X>dam_break_location)+perturbation
    state.q[x_momentum,:,:] = HL*UL*(X<=dam_break_location) + HR*UR*(X>dam_break_location)
    state.q[y_momentum,:,:] = 0.


def step_inclined_source(solver,state,dt):
    """
    source terms for inclined shallow water equation.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    c=(0.7**(3/2) - 1.0)/(0.7 - 1.0)
    dt2=dt/2.0
    q=state.q
    h=q[depth,:,:]
    u=q[x_momentum,:,:]/h+c
    v=q[y_momentum,:,:]/h
    qstar = np.empty(q.shape)
    qstar[0,:,:] = q[0,:,:]
    qstar[1,:,:] = q[1,:,:]+dt2*(h-(u**2+v**2)**(1/2)*u)
    qstar[2,:,:] = q[2,:,:]-dt2*(u**2+v**2)**(1/2)*v
    h=qstar[0,:,:]
    u=qstar[1,:,:]/h+c
    v=qstar[2,:,:]/h

    q[1,:,:] = q[1,:,:]+dt*(h-(u**2+v**2)**(1/2)*u)
    q[2,:,:] = q[2,:,:]-dt*(u**2+v**2)**(1/2)*v

    
def setup(F=2.13,kernel_language='Python', use_petsc=False, outdir='dam_break_F_equals_2_point_13',
          solver_type='classic', riemann_solver='hlle',disable_output=False):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_2D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_2D

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(rs)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1
        solver.step_source = step_inclined_source
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rs)
        solver.step_source = step_inclined_source

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall

    # Domain:
    xlower = 0
    xupper = 20
    mx =2000
    ylower = 0
    yupper = 5
    my = 500
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0/F**2

    qinit(state)

    claw = pyclaw.Controller()
    claw.tfinal = 700
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = 7000
    claw.keep_copy = False

    return claw




if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
