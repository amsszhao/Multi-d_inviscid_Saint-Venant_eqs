function out=extract_sub(f,hr,error_taylor,Rel_error_ODE,Abs_error_ODE,firstcolumn_func,eigen_mat,eigen_vector,eta,r,R,m1,m2,m3)
   t1=datetime('now');
   point=point_sub(f,hr,error_taylor,Rel_error_ODE,Abs_error_ODE,firstcolumn_func,eigen_mat,eigen_vector);
   point=point.classification(eta,r,R,m1,m2,m3);
   out.parameter=[f,hr];
   out.stability=point.stability;
   out.problem=point.problem;
   out.maxchange=point.maxchange;
   t2=datetime('now');
   out.computing_time=t2-t1;
end