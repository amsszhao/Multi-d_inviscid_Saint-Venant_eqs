function out=interior_ode(H,y,ode,F,HR,eigen,eta,lambda)
      out=ode(F,H,HR,eigen,eta,lambda)*y;
end