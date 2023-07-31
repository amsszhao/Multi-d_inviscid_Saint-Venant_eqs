%The following class can be used to study the spectrum of roll waves with
%Froude number F and minimum fluid height Hm

classdef point_roll
    properties
        F % Froude number
        Hm % minimum fluid height
        Hp % maximum fluid height
        error_taylor % control of truncation error of the power series expansion near the sonic point
        ode_rel % relative error used for ode solver
        ode_abs % absolute error used for ode solver
        L % matrices L^{(n)}
        R1 % matrices R^{(1,n)}
        R2 % matrices R^{(2,n)}
        R3 % matrices R^{(3,n)}
        A_func
        First_column %First column of the periodic Evans-Lopatinsky matrix
        maxchange % maximum angle change of contour
        maxchange_allowed % control of maxchange
        stability % stability of roll waves
        low_coefficients % coefficients of low-frequency expansion
        bc % coefficients of low-frequency Weierstrass polynomail 
        low_index %low frequency indexes
        mid_index %low frequency indexes
        problem %detailed stability info
    end
    methods

        % The following function creates a object "roll waves" with Froude
        % number F and minimum fluid height Hm and prepare static matrices
        % for later computations.
        function obj=point_roll(F,Hm,error_taylor,ode_rel,ode_abs,maxchange_allowed)
            obj.F=F;
            obj.Hm=Hm;
            obj.error_taylor=error_taylor;
            obj.ode_rel=ode_rel;
            obj.ode_abs=ode_abs;
            obj.maxchange_allowed=maxchange_allowed;
            obj.Hp=-Hm/2+(Hm^2/4+2/Hm)^0.5;
            % call compute_static_mat to populate matrices L^{(n)},R^{(1,n)},R^{(2,n)}, and R^{(3,n)}
            [obj.L,obj.R1,obj.R2,obj.R3]=compute_static_mat(F);
            % create ode function for evolving interior ode away from the
            % sonic point
            obj.A_func=@(H)[-1/F-1,1,0;H/F^2-(1/F-H*(1/F+1))^2/H^2,-1/F-(2*(1/F-H*(1/F+1)))/H-1,0;0,0,-1/F-(1/F-H*(1/F+1))/H-1];
            obj.maxchange=0;
            % create function that evaluate the first column of the
            % periodic Evans-Lopatinsky determinant
            obj.First_column=@(lam,eta)[lam*(Hm-obj.Hp);obj.Hp-Hm+(lam*(Hm+F*Hm-1))/F-(lam*(obj.Hp+F*obj.Hp-1))/F+(((Hm+F*Hm-1)^2/Hm^2)^(1/2)*(Hm+F*Hm-1))/(F^2*Hm)-(((obj.Hp+F*obj.Hp-1)^2/obj.Hp^2)^(1/2)*(obj.Hp+F*obj.Hp-1))/(F^2*obj.Hp);(eta*(Hm^2*1i-obj.Hp^2*1i))/(2*F^2)];
        end

        % The following function computes the periodic Evan-Lopatinsky at
        % (lam,eta,theta=xi*X)
        function out=periodic_evans_lop(obj,lam,eta,theta)
            n=51;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            % compute W^{(n)}
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);

            % near the sonic point solutions are defined by power series
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            % away from the sonic point solutions are continued by evolving
            % the interior odes
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end

            %compute the periodic Evans-Lopatinsky determinant
            out=det([obj.First_column(lam,eta) Vr1-exp(1i*theta)*Vl1 Vr2-exp(1i*theta)*Vl2]);
        end
        
        % The following function can be used to track the cross-section of
        % the spectra with Re lambda=0 plane
        function out=trace_curve1(obj,x,eta)
            lam=x(1)*1i;
            theta=x(2);
            lop=periodic_evans_lop(obj,lam,eta,theta);
            out=[real(lop),imag(lop)];
        end
        % The following function can be used to track the cross-section of
        % the spectra with Re lambda=0 plane
        function out=trace_curve2(obj,x,theta)
            lam=x(1)*1i;
            eta=x(2);
            lop=periodic_evans_lop(obj,lam,eta,theta);
            out=[real(lop),imag(lop)];
        end

        % The following function can be used to find eta_* when it is
        % called in fminsearch
        function eta=find_eta(obj,theta,eta_guess,lam_i_guess,options)
            out=fsolve(@(x)obj.trace_curve2(x,theta),[lam_i_guess eta_guess],options);
            eta=-out(2);
        end

        % The following function can be used to find local maximum of Re lambda(eta,theta) when it is
        % called in fminsearch
        function lam_real=find_lam_real(obj,x,lam_guess,options)
            eta=x(1);
            theta=x(2);
            lambda=fsolve(@(lambda)obj.search_root_array(lambda,eta,theta),lam_guess,options);
            lam_real=-lambda(1);
        end

        % The following function can be used to find eta_* when a rough
        % (eta_start,lam_i_start,theta_start) is known
        function [eta_a,theta_a,lam_i_a]=find_eta_star(obj,eta_start,lam_i_start,theta_start)
            eta_a=[];
            theta_a=[];
            lam_i_a=[];
            options=optimoptions('fsolve','Algorithm','levenberg-marquardt','FunctionTolerance',1e-7,'StepTolerance',1e-7,'Display','off','UseParallel',true);
            temp_old=fsolve(@(x)obj.trace_curve2(x,theta_start),[lam_i_start eta_start],options);
            eta_a(end+1,:)=temp_old(2);
            theta_a(end+1)=theta_start;
            lam_i_a(end+1)=temp_old(1);
            d_theta=0.02;
            theta_old=theta_start;
            while abs(d_theta)>10^-4
                theta_new=theta_old+d_theta;
                temp_new=fsolve(@(x)obj.trace_curve2(x,theta_new),temp_old,options);
                eta_a(end+1,:)=temp_new(2);
                theta_a(end+1)=theta_new;
                lam_i_a(end+1)=temp_new(1);
                if eta_a(end)<eta_a(end-1)
                    d_theta=d_theta*-1/2;
                end
                theta_old=theta_new;
                temp_old=temp_new;
            end
        end
       
        % The following function can be used to find a spectrum
        % lambda(eta,theta) when it is called by fsolve
        function out=search_root_array(obj,lambda,eta,theta)
            lam=lambda(1)+1i*lambda(2);
            lop=periodic_evans_lop(obj,lam,eta,theta);
            out=[real(lop),imag(lop)];
        end
        
        % The following function computes columns in the periodic
        % Evans-Lopatinsky matrix for fixed lambda and eta
        function [first,Vl1,Vl2,Vr1,Vr2]=get_info_periodic_lop_evans(obj,lam,eta)
            n=51;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end
            first=obj.First_column(lam,eta);
        end
        
        % The following function compute the winding numbers of contours
        % and check stability
        function obj=parallel_verification(obj,path,eta)
            [~,n]=size(path);
            mappath=zeros(1,n);
            first=zeros(3,n);
            Vl1=zeros(3,n);
            Vl2=zeros(3,n);
            Vr1=zeros(3,n);
            Vr2=zeros(3,n);
            parfor i=1:n
                [first(:,i),Vl1(:,i),Vl2(:,i),Vr1(:,i),Vr2(:,i)]=get_info_periodic_lop_evans(obj,path(1,i),eta);
            end
            num_root=0;
            theta_a=0:0.1:2*pi;
            m=length(theta_a);
            j=1;
            while num_root==0 && j<=m
                for i=1:length(path)
                    mappath(i)=det([first(:,i) Vr1(:,i)-exp(1i*theta_a(j))*Vl1(:,i) Vr2(:,i)-exp(1i*theta_a(j))*Vl2(:,i)]);
                end

                % reduce the maximum angle change if it is above allowed
                % value
                [num_root,path,mappath,max_change,first,Vl1,Vl2,Vr1,Vr2]=reduced_maxangle(obj,path,eta,theta_a(j),mappath,first,Vl1,Vl2,Vr1,Vr2);
                if num_root>=1
                    obj.stability=0;
                    obj.problem.num_root=num_root;
                    obj.problem.eta=eta;
                    obj.problem.theta=theta_a(j);
                    obj.problem.mappath=mappath;
                    obj.problem.path=path;
                end
                if obj.maxchange<max_change
                    obj.maxchange=max_change;
                end
                j=j+1;
            end
        end
        % The following function computes additional columns of the
        % periodic Evans-Lopatinsky function at frequencies where the angle change is above allowed value. 
        function [num_root,path,mappath,max_change,first,Vl1,Vl2,Vr1,Vr2]=reduced_maxangle(obj,path,eta,theta,mappath,first,Vl1,Vl2,Vr1,Vr2)
            [num_root,angle_change]=windingnumber(obj,mappath);
            n=sum(abs(angle_change)>obj.maxchange_allowed);
            while n>0
                index_add=zeros(1,n);
                path_add=zeros(2,n);
                first_add=zeros(3,n);
                Vl1_add=zeros(3,n);
                Vl2_add=zeros(3,n);
                Vr1_add=zeros(3,n);
                Vr2_add=zeros(3,n);
                mappath_add=zeros(1,n);
                k=1;
                for i=1:length(angle_change)
                    if abs(angle_change(i))>obj.maxchange_allowed
                        index_add(k)=i;
                        path_add(2,k)=path(2,i);
                        if path_add(2,k)==0
                            path_add(1,k)=(path(1,i)+path(1,i+1))/2;
                        else
                            path_add(1,k)=abs(path(1,i))*exp(1i*(angle(path(1,i))+angle(path(1,1+i)))/2);
                        end
                        k=k+1;
                    end
                end
                parfor i=1:n
                    [first_add(:,i),Vl1_add(:,i),Vl2_add(:,i),Vr1_add(:,i),Vr2_add(:,i)]=get_info_periodic_lop_evans(obj,path_add(1,i),eta);
                    mappath_add(i)=det([first_add(:,i) Vr1_add(:,i)-exp(1i*theta)*Vl1_add(:,i) Vr2_add(:,i)-exp(1i*theta)*Vl2_add(:,i)]);
                end
                for i=1:n
                    first=[first(:,1:(index_add(i)+i-1)) first_add(:,i) first(:,(index_add(i)+i):end)];
                    Vl1=[Vl1(:,1:(index_add(i)+i-1)) Vl1_add(:,i) Vl1(:,(index_add(i)+i):end)];
                    Vl2=[Vl2(:,1:(index_add(i)+i-1)) Vl2_add(:,i) Vl2(:,(index_add(i)+i):end)];
                    Vr1=[Vr1(:,1:(index_add(i)+i-1)) Vr1_add(:,i) Vr1(:,(index_add(i)+i):end)];
                    Vr2=[Vr2(:,1:(index_add(i)+i-1)) Vr2_add(:,i) Vr2(:,(index_add(i)+i):end)];
                    path=[path(:,1:(index_add(i)+i-1)) path_add(:,i) path(:,(index_add(i)+i):end)];
                    mappath=[mappath(1:(index_add(i)+i-1)) mappath_add(:,i) mappath(:,(index_add(i)+i):end)];
                end
                [num_root,angle_change]=windingnumber(obj,mappath);
                n=sum(abs(angle_change)>obj.maxchange_allowed);
            end
            max_change=max(abs(angle_change));
            
        end

        % The following function computes partial_theta Delta(lambda,eta,0)
        function out=periodic_lop_evans1(obj,lam,eta)
            n=21;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end
            first=obj.First_column(lam,eta);
            out=Vl1(1)*Vl2(2)*first(3)*2i - Vl1(1)*Vl2(3)*first(2)*2i - Vl1(2)*Vl2(1)*first(3)*2i + Vl1(2)*Vl2(3)*first(1)*2i + Vl1(3)*Vl2(1)*first(2)*2i - Vl1(3)*Vl2(2)*first(1)*2i - Vl1(1)*Vr2(2)*first(3)*1i + Vl1(1)*Vr2(3)*first(2)*1i + Vl1(2)*Vr2(1)*first(3)*1i - Vl1(2)*Vr2(3)*first(1)*1i - Vl1(3)*Vr2(1)*first(2)*1i + Vl1(3)*Vr2(2)*first(1)*1i + Vl2(1)*Vr1(2)*first(3)*1i - Vl2(1)*Vr1(3)*first(2)*1i - Vl2(2)*Vr1(1)*first(3)*1i + Vl2(2)*Vr1(3)*first(1)*1i + Vl2(3)*Vr1(1)*first(2)*1i - Vl2(3)*Vr1(2)*first(1)*1i;
        end

        % The following function computes partial_theta^2 Delta(lambda,eta,0)
        function out=periodic_lop_evans2(obj,lam,eta)
            n=21;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end
            first=obj.First_column(lam,eta);
            out=4*Vl1(1)*Vl2(3)*first(2) - 4*Vl1(1)*Vl2(2)*first(3) + 4*Vl1(2)*Vl2(1)*first(3) - 4*Vl1(2)*Vl2(3)*first(1) - 4*Vl1(3)*Vl2(1)*first(2) + 4*Vl1(3)*Vl2(2)*first(1) + Vl1(1)*Vr2(2)*first(3) - Vl1(1)*Vr2(3)*first(2) - Vl1(2)*Vr2(1)*first(3) + Vl1(2)*Vr2(3)*first(1) + Vl1(3)*Vr2(1)*first(2) - Vl1(3)*Vr2(2)*first(1) - Vl2(1)*Vr1(2)*first(3) + Vl2(1)*Vr1(3)*first(2) + Vl2(2)*Vr1(1)*first(3) - Vl2(2)*Vr1(3)*first(1) - Vl2(3)*Vr1(1)*first(2) + Vl2(3)*Vr1(2)*first(1);
        end

        % The following function computes partial_theta^3 Delta(lambda,eta,0)
        function out=periodic_lop_evans3(obj,lam,eta)
            n=21;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end
            first=obj.First_column(lam,eta);
            out=- Vl1(1)*Vl2(2)*first(3)*8i + Vl1(1)*Vl2(3)*first(2)*8i + Vl1(2)*Vl2(1)*first(3)*8i - Vl1(2)*Vl2(3)*first(1)*8i - Vl1(3)*Vl2(1)*first(2)*8i + Vl1(3)*Vl2(2)*first(1)*8i + Vl1(1)*Vr2(2)*first(3)*1i - Vl1(1)*Vr2(3)*first(2)*1i - Vl1(2)*Vr2(1)*first(3)*1i + Vl1(2)*Vr2(3)*first(1)*1i + Vl1(3)*Vr2(1)*first(2)*1i - Vl1(3)*Vr2(2)*first(1)*1i - Vl2(1)*Vr1(2)*first(3)*1i + Vl2(1)*Vr1(3)*first(2)*1i + Vl2(2)*Vr1(1)*first(3)*1i - Vl2(2)*Vr1(3)*first(1)*1i - Vl2(3)*Vr1(1)*first(2)*1i + Vl2(3)*Vr1(2)*first(1)*1i;
        end
        
        % The following function computes partial_theta^4 Delta(lambda,eta,0)
        function out=periodic_lop_evans4(obj,lam,eta)
            n=21;
            optionODE=odeset('RelTol',obj.ode_rel,'AbsTOl',obj.ode_abs);
            W_con=compute_W(n,lam,eta,obj.F,obj.L,obj.R1,obj.R2,obj.R3);
            temp=abs(W_con(:,:,n)./W_con(:,:,1));
            max_temp1=0;
            if temp(1,1)~=Inf && temp(1,1)>max_temp1
                max_temp1=temp(1,1);
            end
            if temp(2,1)~=Inf && temp(2,1)>max_temp1
                max_temp1=temp(2,1);
            end
            if temp(3,1)~=Inf && temp(3,1)>max_temp1
                max_temp1=temp(3,1);
            end
            max_temp2=0;
            if temp(1,2)~=Inf && temp(1,2)>max_temp2
                max_temp2=temp(1,2);
            end
            if temp(2,2)~=Inf && temp(2,2)>max_temp2
                max_temp2=temp(2,2);
            end
            if temp(3,2)~=Inf && temp(3,2)>max_temp2
                max_temp2=temp(3,2);
            end
            dh1=(obj.error_taylor/max_temp1)^(1/(n-1));
            dh2=(obj.error_taylor/max_temp2)^(1/(n-1));
            Wr1=W_con(:,1,1);
            if dh1>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh^(i-1);
                end
                Vr1=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
            else
                for i=2:n
                    Wr1=Wr1+W_con(:,1,i)*dh1^(i-1);
                end
                Vr1=obj.A_func(1+dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh1 obj.Hp],Vr1,optionODE);
                Vr1=transpose(out(end,:));
            end
            
            Wr2=W_con(:,2,1);
            if dh2>(obj.Hp-1)
                dh=(obj.Hp-1);
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh^(i-1);
                end
                Vr2=obj.A_func(obj.Hp)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
            else
                for i=2:n
                    Wr2=Wr2+W_con(:,2,i)*dh2^(i-1);
                end
                Vr2=obj.A_func(1+dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wr2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1+dh2 obj.Hp],Vr2,optionODE);
                Vr2=transpose(out(end,:));
            end
            Wl1=W_con(:,1,1);
            if (1-dh1)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*dh^(i-1);
                end
                Vl1=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
            else
                for i=2:n
                    Wl1=Wl1+W_con(:,1,i)*(-dh1)^(i-1);
                end
                Vl1=obj.A_func(1-dh1)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl1;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh1 obj.Hm],Vl1,optionODE);
                Vl1=transpose(out(end,:));
            end
            
            Wl2=W_con(:,2,1);
            if (1-dh2)<obj.Hm
                dh=(obj.Hm-1);
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*dh^(i-1);
                end
                Vl2=obj.A_func(obj.Hm)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
            else
                for i=2:n
                    Wl2=Wl2+W_con(:,2,i)*(-dh2)^(i-1);
                end
                Vl2=obj.A_func(1-dh2)*[1,0,obj.F/(obj.F+1);0,0,1;0,1,0]*Wl2;
                [~,out]=ode23tb(@(H,y)interior_ode(H,y,lam,eta,obj.F),[1-dh2 obj.Hm],Vl2,optionODE);
                Vl2=transpose(out(end,:));
            end
            first=obj.First_column(lam,eta);
            out=16*Vl1(1)*Vl2(2)*first(3) - 16*Vl1(1)*Vl2(3)*first(2) - 16*Vl1(2)*Vl2(1)*first(3) + 16*Vl1(2)*Vl2(3)*first(1) + 16*Vl1(3)*Vl2(1)*first(2) - 16*Vl1(3)*Vl2(2)*first(1) - Vl1(1)*Vr2(2)*first(3) + Vl1(1)*Vr2(3)*first(2) + Vl1(2)*Vr2(1)*first(3) - Vl1(2)*Vr2(3)*first(1) - Vl1(3)*Vr2(1)*first(2) + Vl1(3)*Vr2(2)*first(1) + Vl2(1)*Vr1(2)*first(3) - Vl2(1)*Vr1(3)*first(2) - Vl2(2)*Vr1(1)*first(3) + Vl2(2)*Vr1(3)*first(1) + Vl2(3)*Vr1(1)*first(2) - Vl2(3)*Vr1(2)*first(1);
        end
        
        % The following function computes windingnumber of the contour
        % mappath
        function [num_root,angle_change]=windingnumber(~,mappath)
            nn=length(mappath)-1;
            real_mappath=real(mappath);
            imag_mappath=imag(mappath);
            angle_change=sign(real_mappath(1:nn).*imag_mappath(2:nn+1)-imag_mappath(1:nn).*real_mappath(2:nn+1)).*acos((real_mappath(1:nn).*real_mappath(2:nn+1)+imag_mappath(1:nn).*imag_mappath(2:nn+1))./abs(mappath(1:nn))./abs(mappath(2:nn+1)));
            num_root=sum(angle_change);
            num_root=round(num_root/2/pi);
        end

        %The following function compute the integrand of low frequency
        %expansion coefficients a_{(k1,k2,0)}
        function out=integrand0(obj,theta1,theta2,k1,k2,r)
            [n1,n2]=size(theta1);
            out=zeros([n1,n2]);
            parfor i=1:n1*n2
                out(i)=periodic_evans_lop(obj,r*exp(1i*theta1(i)),r*exp(1i*theta2(i)),0)/(4*pi^2*r^(k1+k2)*exp(1i*theta1(i)*k1+1i*theta2(i)*k2));
            end
        end

        %The following function compute the integrand of low frequency
        %expansion coefficients a_{(k1,k2,1)}
        function out=integrand1(obj,theta1,theta2,k1,k2,r)
            [n1,n2]=size(theta1);
            out=zeros([n1,n2]);
            parfor i=1:n1*n2
                out(i)=periodic_lop_evans1(obj,r*exp(1i*theta1(i)),r*exp(1i*theta2(i)))/(4*pi^2*r^(k1+k2)*exp(1i*theta1(i)*k1+1i*theta2(i)*k2));
            end
        end

        %The following function compute the integrand of low frequency
        %expansion coefficients a_{(k1,k2,2)}
        function out=integrand2(obj,theta1,theta2,k1,k2,r)
            [n1,n2]=size(theta1);
            out=zeros([n1,n2]);
            parfor i=1:n1*n2
                out(i)=periodic_lop_evans2(obj,r*exp(1i*theta1(i)),r*exp(1i*theta2(i)))/(2*4*pi^2*r^(k1+k2)*exp(1i*theta1(i)*k1+1i*theta2(i)*k2));
            end
        end

        %The following function compute the integrand of low frequency
        %expansion coefficients a_{(k1,k2,3)}
        function out=integrand3(obj,theta1,theta2,k1,k2,r)
            [n1,n2]=size(theta1);
            out=zeros([n1,n2]);
            parfor i=1:n1*n2
                out(i)=periodic_lop_evans3(obj,r*exp(1i*theta1(i)),r*exp(1i*theta2(i)))/(6*4*pi^2*r^(k1+k2)*exp(1i*theta1(i)*k1+1i*theta2(i)*k2));
            end
        end

        %The following function compute the integrand of low frequency
        %expansion coefficients a_{(k1,k2,4)}
        function out=integrand4(obj,theta1,theta2,k1,k2,r)
            [n1,n2]=size(theta1);
            out=zeros([n1,n2]);
            parfor i=1:n1*n2
                out(i)=periodic_lop_evans4(obj,r*exp(1i*theta1(i)),r*exp(1i*theta2(i)))/(24*4*pi^2*r^(k1+k2)*exp(1i*theta1(i)*k1+1i*theta2(i)*k2));
            end
        end

        %The following function compute the integrand of mid frequency
        %expansion coefficients a_{(1,0,0)}
        function out=mid100(obj,lamstar,thetastar,theta,r)
            n=length(theta);
            out=zeros(1,n);
            parfor i=1:n
                out(i)=periodic_evans_lop(obj,lamstar+r*exp(1i*theta(i)),0,thetastar)/(2*pi*r*exp(1i*theta(i)));
            end
        end

        %The following function compute the integrand of mid frequency
        %expansion coefficients a_{(0,2,0)}
        function out=mid020(obj,lamstar,thetastar,theta,r)
            n=length(theta);
            out=zeros(1,n);
            parfor i=1:n
                out(i)=periodic_evans_lop(obj,lamstar,r*exp(1i*theta(i)),thetastar)/(2*pi*r^2*exp(2i*theta(i)));
            end
        end

        %The following function generates a half_annulus in the lambda
        %frequency plane centered at the origin
        function path=setcontour_half_annulus(~,r,R,m1,m2,m3)
            path1=r*exp(1i*(pi/2:-pi/m1:-pi/2));
            path2=(-r:(r-R)/m2:-R)*1i;
            path3=R*exp(1i*(-pi/2:pi/m3:pi/2));
            path4=(R:(r-R)/m2:r)*1i;
            path=[[path1,path2(2:m2+1),path3(2:m3+1),path4(2:m2+1)];[ones(1,m1) zeros(1,m2) ones(1,m3) zeros(1,1+m2)]];
        end
        
        %The following function generates a circle in the lambda
        %frequency plane 
        function path=setcontour_circle(~,center,r,n)
            path=center+exp(1i*(0:2*pi/n:2*pi))*r;
        end
        
        %The following function generates a half circle in the lambda
        %frequency plane centered at the origin
        function path=setcontour_half_circle(~,r,R,m1,m2,m3)
            path1=(r:-(2*r)/m1:(-r))*1i;
            path2=(-r:(r-R)/m2:-R)*1i;
            path3=R*exp(1i*(-pi/2:pi/m3:pi/2));
            path4=(R:(r-R)/m2:r)*1i;
            path=[[path1,path2(2:m2+1),path3(2:m3+1),path4(2:m2+1)];[zeros(1,m1+m2) ones(1,m3) zeros(1,1+m2)]];
        end

        % The following function computes low frequency expansion
        % coefficents, coeffients of the Weierstrass polynomial and
        % low-frequency indexes
        function obj=low_frequency(obj)
            a120=integral2(@(theta1,theta2)obj.integrand0(theta1,theta2,1,2,0.01),0,2*pi,0,2*pi);
            a300=integral2(@(theta1,theta2)obj.integrand0(theta1,theta2,3,0,0.01),0,2*pi,0,2*pi);
            a101=integral2(@(theta1,theta2)obj.integrand1(theta1,theta2,1,0,0.01),0,2*pi,0,2*pi);
            a020=integral2(@(theta1,theta2)obj.integrand0(theta1,theta2,0,2,0.01),0,2*pi,0,2*pi);
            a200=integral2(@(theta1,theta2)obj.integrand0(theta1,theta2,2,0,0.01),0,2*pi,0,2*pi);
            a021=integral2(@(theta1,theta2)obj.integrand1(theta1,theta2,0,2,0.01),0,2*pi,0,2*pi);
            a201=integral2(@(theta1,theta2)obj.integrand1(theta1,theta2,2,0,0.01),0,2*pi,0,2*pi);
            a102=integral2(@(theta1,theta2)obj.integrand2(theta1,theta2,1,0,0.01),0,2*pi,0,2*pi);
            b01=imag(a101/a200)*1i;
            c20=real(a020/a200);
            c21=imag((a021 - (a020*(a201 - (a101*a300)/a200))/a200)/a200)*1i;
            b20=real(-(a020*a300 - a120*a200)/a200^2);
            b02=real((a300*a101^2 - a201*a101*a200 + a102*a200^2)/a200^3);
            obj.bc=[b01 b02 b20 c20 c21];
            obj.low_coefficients=[a200 a101 a020 a120 a300 a021 a201 a102];
            obj.low_index=[b02 b20 c21^2*(b01^2*b20^2-2*b01*b20*c21+4*b02*c20*b20+c21^2) c21^2-b01*b20*c21+2*b02*b20*c20 c20*b02-b01*c21];
        end


        % The following function computes mid frequency expansion
        % coefficents, coeffients of the Weierstrass polynomial and
        % mid-frequency indexes
        function obj=mid_frequency(obj,lamstar,thetastar)
            a100=2*integral(@(theta)obj.mid100(lamstar,thetastar,theta,min([0.1,(obj.F-2)/5,imag(lamstar)/4])),0,pi,'AbsTol',1e-5,'RelTol',1e-5);
            a020=2*integral(@(theta)obj.mid020(lamstar,thetastar,theta,min([0.1,(obj.F-2)/5,imag(lamstar)/4])),0,pi,'AbsTol',1e-5,'RelTol',1e-5);
            obj.mid_index=[-a020/a100 a100 a020];
        end
        
        % The following function examines mid frequency stability 
        function obj=classification(obj,eta,r,R,m1,m2,m3)
            obj.stability=1;
            mm=length(r);
            nn=length(eta);
            i=1;
            while (i<=mm)&&obj.stability
                path=setcontour_half_annulus(obj,r(i),R(i),m1(i),m2(i),m3(i));
                obj=parallel_verification(obj,path,0);
                i=i+1;
            end
            for j=2:nn
                i=1;
                while (i<=mm)&&obj.stability
                    path=setcontour_half_circle(obj,r(i),R(i),m1(i),m2(i),m3(i));
                    obj=parallel_verification(obj,path,eta(j));
                    i=i+1;
                end
            end
        end
    end
end