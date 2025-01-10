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

def qinit(state,f,xin,period,num_mesh,num_period):
    from scipy.integrate import odeint
    fun=lambda h,t:(f**2*h**2-(1.0+2.0*f)*h+1.0)/(h**2+h+1.0)
    t=np.linspace(xin[0],xin[0]+period[0],num_mesh)
    tn=t[t<0]
    tp=t[t>0]
    tn=np.append([0],tn[::-1])
    tp=np.append([0],tp)
    sol1=odeint(fun,1,tn)
    sol2=odeint(fun,1,tp)
    sol1=sol1[1:]
    sol2=sol2[1:]
    sol=np.append(sol1[::-1],sol2)
    harray=np.tile(sol,num_period)
    huarray=-1.0/f
    
    X, Y = state.p_centers
    perturbation=np.zeros(X.shape)
    stateq=np.zeros(state.q.shape)
    
    for j in range(100):
        stateq[0,:,j] = harray
        stateq[1,:,j] = huarray
        stateq[2,:,j] = 0.0
    	    
    for i in range(8000):
        for j in range(100):
            r2=(X[i,j]-5)**2 + (Y[i,j]-0.08)**2*2
            if r2<0.01:
                perturbation[i,j]=np.exp(-0.03/(0.01-r2))
    state.q[depth     ,:,:]=stateq[0,:,:]+perturbation
    state.q[x_momentum,:,:] = stateq[1,:,:]
    state.q[y_momentum,:,:] = 0.0


def step_inclined_source(solver,state,dt):
    """
    source terms for inclined shallow water equation.
    Integrated using a 2-stage, 2nd-order Runge-Kutta method.
    This is a Clawpack-style source term routine, which approximates
    the integral of the source terms over a step.
    """
    c=1.0+1.0/6.0
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
    
    
def setup(f=6.0,hn=0.28,num_mesh=800,num_period=10,kernel_language='Python', use_petsc=False, outdir='roll_width_point16',
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
        solver.limiters = pyclaw.limiters.tvd.vanleer
        solver.dimensional_split=1
        solver.step_source = step_inclined_source
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rs)
        solver.step_source = step_inclined_source

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    
    from scipy import integrate
    fun=lambda h:(h**2+h+1.0)/(f**2*h**2-(1.0+2.0*f)*h+1.0)
    period=integrate.quad(fun,hn,-hn/2.0+(hn**2/4.0+2.0/hn)**0.5)
    xin=integrate.quad(fun,1,hn)
    xlower = 0.0
    xupper = num_period*period[0]
    mx =num_period*num_mesh
    ylower = 0
    yupper = 0.16
    my = 100
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0/f**2

    qinit(state,f,xin,period,num_mesh,num_period)

    claw = pyclaw.Controller()
    claw.tfinal = 1000
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = 5000
    claw.keep_copy = False

    return claw




if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup)
