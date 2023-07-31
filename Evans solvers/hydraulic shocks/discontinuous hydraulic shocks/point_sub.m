%The following class can be used to study the spectrum of discontinuous hydraulic shocks with
%Froude number F and right fluid height HR
classdef point_sub
    properties
        F % Froude number
        HR % right fluid height
        Hs % fluid height after the subshock H_*
        error_taylor % control of truncation error of the power series expansion near -infty
        Rel_error_ODE % relative error used for ode solver
        Abs_error_ODE % absolute error used for ode solver
        L % matrices L_i
        R1 % matrices R_{1i}
        R2 % matrices R_{2i}
        R3 % matrices R_{3i}
        R4 % matrices R_{4i}
        firstcolumn_func %First column of the Evans-Lopatinsky matrix
        eigen_mat % eigenvalue of limiting matrix
        eigen_vector % eigenvector of limiting matrix
        maxchange % max angle change
        stability % stability of discontinuous hydraulic shocks
        problem % detail of stability info
    end
    methods
        % The following function creates a object "a discontinuous hydraulic shock" with Froude
        % number F and right fluid height HR and prepare static matrices
        % for later computations.
        function obj=point_sub(f,hr,error_taylor,Rel_error_ODE,Abs_error_ODE,firstcolumn_func,eigen_mat,eigen_vector)
            obj.F=f;
            obj.HR=hr;
            obj.Hs=-(hr - (hr*(hr + 8*f^2 + hr^2 + 2*hr^(3/2)))^(1/2) + hr^(3/2))/(2*hr^(1/2) + 2);
            obj.error_taylor=error_taylor;
            obj.Rel_error_ODE=Rel_error_ODE;
            obj.Abs_error_ODE=Abs_error_ODE;
            obj.L=zeros(6,1);
            obj.R1=cell(6,1);
            obj.R2=cell(7,1);
            obj.R3=cell(7,1);
            obj.R4=zeros(6,1);
            % call script file static_constant_v to compute for static
            % matrices
            static_constant_v
            % create function that evaluate the first column of the
            % Evans-Lopatinsky determinant
            obj.firstcolumn_func=firstcolumn_func;
            obj.eigen_mat=eigen_mat;
            obj.eigen_vector=eigen_vector;
            obj.maxchange=0;
        end

        % The following function computes the a non-radial Evans-Lopatinsky at
        % (lam,eta=0)
        function out=lopatinsky0_non_radial(obj,lam)
            n=50;
            compute_V0_n
            temp=abs(V(:,n)./V(:,1));
            max_temp=0;
            if temp(1)~=Inf && temp(1)>max_temp
                max_temp=temp(1);
            end
            if temp(2)~=Inf && temp(2)>max_temp
                max_temp=temp(2);
            end
            if temp(3)~=Inf && temp(3)>max_temp
                max_temp=temp(3);
            end
            dh=sign(obj.Hs-1)*min((obj.error_taylor/max_temp)^(1/(n-1)),abs(1-obj.Hs));
            V_L=V(:,1);
            for i=2:n
                V_L=V_L+V(:,i)*dh^(i-1);
            end
            V_L=V_L/norm(V_L);
            if (1+dh)~=obj.Hs
                optionODE=odeset('RelTol',obj.Rel_error_ODE,'AbsTOl',obj.Abs_error_ODE);
                [~,yy]=ode23tb(@(H,y)interior0_ode_nonradial(H,y,obj.F,obj.HR,lam),[1+dh obj.Hs],V_L,optionODE);
                V_L=transpose(yy(end,:));
            end
            first_col=obj.firstcolumn_func(obj.F,obj.HR,0,lam);
            out=V_L(1)*first_col(1)+V_L(2)*first_col(2)+V_L(3)*first_col(3);
        end
        
        % The following function computes the a non-radial Evans-Lopatinsky at
        % (lam,eta ~= 0)
        function out=lopatinsky_non_radial(obj,lam,eta)
            n=50;
            compute_V_n
            temp=abs(V(:,n)./V(:,1));
            max_temp=0;
            if temp(1)~=Inf && temp(1)>max_temp
                max_temp=temp(1);
            end
            if temp(2)~=Inf && temp(2)>max_temp
                max_temp=temp(2);
            end
            if temp(3)~=Inf && temp(3)>max_temp
                max_temp=temp(3);
            end
            dh=sign(obj.Hs-1)*min((obj.error_taylor/max_temp)^(1/(n-1)),abs(1-obj.Hs));
            V_L=V(:,1);
            for i=2:n
                V_L=V_L+V(:,i)*dh^(i-1);
            end
            V_L=V_L/norm(V_L);
            if (1+dh)~=obj.Hs
                optionODE=odeset('RelTol',obj.Rel_error_ODE,'AbsTOl',obj.Abs_error_ODE);
                [~,yy]=ode23tb(@(H,y)interior_ode_nonradial(H,y,obj.F,obj.HR,eta,lam),[1+dh obj.Hs],V_L,optionODE);
                V_L=transpose(yy(end,:));
            end
            first_col=obj.firstcolumn_func(obj.F,obj.HR,eta,lam);
            out=V_L(1)*first_col(1)+V_L(2)*first_col(2)+V_L(3)*first_col(3);
        end
        %The following function generates a half_annulus in the lambda
        %frequency plane centered at the origin
        function path=setcontour_half_annulus(~,r,R,m1,m2,m3)
            path1=r*exp(1i*(pi/2:-pi/m1:-pi/2));
            path2=(-r:(r-R)/m2:-R)*1i;
            path3=R*exp(1i*(-pi/2:pi/m3:pi/2));
            path4=(R:(r-R)/m2:r)*1i;
            path=[path1,path2(2:m2+1),path3(2:m3+1),path4(2:m2+1)];
        end
        %The following function generates a half circle in the lambda
        %frequency plane centered at the origin
        function path=setcontour_half_circle(~,r,R,m1,m2,m3)
            path1=(r:-(2*r)/m1:(-r))*1i;
            path2=(-r:(r-R)/m2:-R)*1i;
            path3=R*exp(1i*(-pi/2:pi/m3:pi/2));
            path4=(R:(r-R)/m2:r)*1i;
            path=[path1,path2(2:m2+1),path3(2:m3+1),path4(2:m2+1)];
        end

        % The function below compute Evans-Lopatinsky with eta=0 on the lambda-path
        function mappath=parallel_lop0(obj,path)
            pathlen=length(path);
            mappath=zeros(pathlen,1);
            parfor i=1:pathlen
                mappath(i)=lopatinsky0_non_radial(obj,path(i));
            end
        end

        % The function below compute Evans-Lopatinsky with eta~=0 on the lambda-path
        function mappath=parallel_lop(obj,path,eta)
            pathlen=length(path);
            mappath=zeros(pathlen,1);
            parfor i=1:pathlen
                mappath(i)=lopatinsky_non_radial(obj,path(i),eta);
            end
        end

        % The following function computes windingnumber of the contour
        % mappath
        function [num_root,max_change]=windingnumber(~,mappath)
            nn=length(mappath)-1;
            real_mappath=real(mappath);
            imag_mappath=imag(mappath);
            angle_change=sign(real_mappath(1:nn).*imag_mappath(2:nn+1)-imag_mappath(1:nn).*real_mappath(2:nn+1)).*acos((real_mappath(1:nn).*real_mappath(2:nn+1)+imag_mappath(1:nn).*imag_mappath(2:nn+1))./abs(mappath(1:nn))./abs(mappath(2:nn+1)));
            num_root=sum(angle_change);
            max_change=max(abs(angle_change));
            num_root=[round(num_root/2/pi) num_root/2/pi];
        end

        % The following function examines mid frequency stability 
        function obj=classification(obj,eta,r,R,m1,m2,m3)
            obj.stability=1;
            mm=length(r);
            nn=length(eta);
            i=1;
            while (i<=mm)&&obj.stability
                path=setcontour_half_annulus(obj,r(i),R(i),m1(i),m2(i),m3(i));
                mappath=parallel_lop0(obj,path);
                [num_root,max_change]=windingnumber(obj,mappath);
                if max_change>obj.maxchange
                    obj.maxchange=max_change;
                end
                if num_root(1)>0
                    obj.stability=false;
                    obj.problem.max_change=max_change;
                    obj.problem.circle=[0 r(i) R(i)];
                    obj.problem.num_root=num_root;
                    obj.problem.mappath=mappath;
                    break;
                end
		    i=i+1;
            end
            for j=2:nn
                i=1;
                while (i<=mm)&&obj.stability
                    path=setcontour_half_circle(obj,r(i),R(i),m1(i),m2(i),m3(i));
                    mappath=parallel_lop(obj,path,eta(j));
                    [num_root,max_change]=windingnumber(obj,mappath);
                    if max_change>obj.maxchange
                        obj.maxchange=max_change;
                    end
                    if num_root(1)>0
                        obj.stability=false;
                        obj.problem.max_change=max_change;
                        obj.problem.circle=[eta(j) r(i) R(i)];
                        obj.problem.num_root=num_root;
                        obj.problem.mappath=mappath;
                        break;
                    end
                    i=i+1;
                end
            end
        end
    end
end
