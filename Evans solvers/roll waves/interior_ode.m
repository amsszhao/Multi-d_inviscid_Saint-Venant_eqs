%The following function computes the interior odes
function out=interior_ode(H,y,lam,eta,F)
A=[-1/F-1,1,0;H/F^2-(1/F-H*(1/F+1))^2/H^2,-1/F-(2*(1/F-H*(1/F+1)))/H-1,0;0,0,-1/F-(1/F-H*(1/F+1))/H-1];
E=[0,0,0;(F^2*H^3+2*F^2*H^2+4*F*H^2-4*F*H+2*H^2-4*H+2)/(F^2*H^3),-(2*(H+F*H-1))/(F*H^2),0;0,0,-(H+F*H-1)/(F*H^2)];
A2=[0,0,1;0,0,(H+F*H-1)/(F*H);H/F^2,0,0];
Hpinverse=-(H^2 + H + 1)/(- F^2*H^2 + 2*F*H + H - 1);
out=Hpinverse*(E-lam*eye(3)-1i*eta*A2)*A^-1*y;
end

