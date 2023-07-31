function out=interior0_ode(H,y,ode,F,HR,eigen,lambda)
      out=ode(F,H,HR,eigen,lambda)*y;
end