function out=interior0_ode_nonradial(H,y,F,HR,lambda)
E=[0,0,0;(2*(H-HR+H*HR^(1/2)+H*HR)^2)/(H^3*(HR^(1/2)+1)^2)+1,-(2*(H-HR+H*HR^(1/2)+H*HR))/(H^2*(HR^(1/2)+1)),0;0,0,-(H-HR+H*HR^(1/2)+H*HR)/(H^2*(HR^(1/2)+1))];
A=[-(HR^(3/2)-1)/(HR-1),1,0;H/F^2-((HR-HR^(3/2))/(HR-1)+(H*(HR^(3/2)-1))/(HR-1))^2/H^2,(H-2*HR+H*HR^(1/2)+H*HR)/(H*(HR^(1/2)+1)),0;0,0,-HR/(H*(HR^(1/2)+1))];
Hprime=(F^2*H^2*(HR^(1/2)+1)^2*(H-(((H-HR+H*HR^(1/2)+H*HR)^2/H^2)^(1/2)*(H-HR+H*HR^(1/2)+H*HR))/(H*(HR^(1/2)+1)^2)))/(H^3*HR+H^3-F^2*HR^2+2*H^3*HR^(1/2));
mat=-transpose((E-lambda*eye(3))*A^-1)/Hprime;
out=mat*y-y*y'*mat*y;
end