% The following function extract stability info for the roll wave with
% Froude number F and minimum fluid height Hm.
function out=extract_roll(F,Hm,error_taylor,ode_rel,ode_abs,maxchange_allowed,eta,r,R,m1,m2,m3)
    t1=datetime('now');
    %create object "roll waves" with Froude number F and minimum fluid
    %height Hm
    point=point_roll(F,Hm,error_taylor,ode_rel,ode_abs,maxchange_allowed);
    
    % examine mid-frequency stability 
    point=point.classification(eta,r,R,m1,m2,m3);

    % compute low-frequency indexes
    point=low_frequency(point);

    % extract stability info for records
    out.parameter=[F,Hm];
    out.stability=point.stability;
    out.problem=point.problem;
    out.maxchange=point.maxchange;
    out.low_index=point.low_index;
    t2=datetime('now');
    out.computing_time=t2-t1;
end
