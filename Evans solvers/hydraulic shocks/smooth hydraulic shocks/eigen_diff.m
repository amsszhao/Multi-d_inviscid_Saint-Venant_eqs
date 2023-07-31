function out=eigen_diff(x,xi,F_v,HR_v,gamma_sol_func,coeff_eq)
gamma_v=gamma_sol_func(F_v,HR_v,x(1)+x(2)*1i,xi);
coeff_v=coeff_eq(F_v,HR_v,x(1)+x(2)*1i,xi);
temp=coeff_v(1)*gamma_v^2+coeff_v(2)*gamma_v+coeff_v(3);
out=[real(temp) imag(temp)];
end

 