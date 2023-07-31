clc
clear
filepath=strcat(pwd,'/roll_width_point18_periodic');
N=5000;
mov(1:N)= struct('cdata',[],'colormap',[]);
filenamev='roll_width_point18_periodic';

hn=0.28;
f=6;
num_mesh=800;
num_period=10;
period=integral(@(h)(h.^2+h+1.0)./(f^2*h.^2-(1.0+2.0*f)*h+1.0),hn,-hn/2.0+(hn^2/4.0+2.0/hn)^0.5);
dx=period/800;
X_lower=0;
X_upper=period*num_period-dx;
dx=period/num_mesh;
Y_lower=0.000;
Y_upper=0.18;
dy=0.18/99;
[X,Y] = meshgrid(X_lower:dx:X_upper,Y_lower:dy:Y_upper);
[m,n]=size(X);

dt=0.2;
f=figure;
for i=1:N
    if i-1<=9
        format='fort.q000%d';
    else if i-1<=99
            format='fort.q00%d';
        else if i-1<=999
                format='fort.q0%d';
            else
                format='fort.q%d';
            end
        end
    end
    filename=fullfile(filepath,sprintf(format,i-1));
    fileID=fopen(filename);
    C=textscan(fileID,'%f %f %f','HeaderLines',9);
    fclose(fileID);
    H=C{1};
    H=reshape(H,[n,m]);
    H=transpose(H);
    shock_2d=zeros(100,10);
    k=1;
    for l=1:100
        kk=1;
        dh=H(l,2:end)-H(l,1:end-1);
        flag=1;
        sum_dh=[];
        shock_index=0;
        for j=2:7999
            if dh(j)>=0
                flag=1;
                sum_dh(end+1)=dh(j);
            else
                if (flag==1) && (dh(j)<-0.5) && (j>(shock_index+5))
                    shock(k,:)=[X(1,j) Y(l,1)];
                    shock_2d(l,kk)=X(1,j);
                    kk=kk+1;
                    k=k+1;
                    flag=-1;
                    sum_dh=[];
                    shock_index=j;
                end
            end
        end
    end
    ax=subplot(2,1,1);
    s=surf(X,Y,H);
    s.EdgeColor='none';
    view(3)
    ax.XLim = [X_lower X_upper];
    ax.YLim = [0 0.18];
    ax.CameraPosition =[(X_lower+X_upper)/2 -0.09 5];
    ax.CameraTarget = [(X_lower+X_upper)/2 0.09 1.5];
    ax.CameraUpVector = [0 0 1];
    ax.CameraViewAngle = 35;
    ax.DataAspectRatio = [0.4 0.4 1];
    l1 = light;
    theta=40/180*pi;
    alpha=5/180*pi;
    l1.Position = [sin(theta)*sin(alpha) sin(theta)*cos(alpha) cos(theta)];
    l1.Style ='infinite';
    l1.Color = [249/256 247/256 248/256];
    s.FaceColor = [0 105/256 148/256];
    s.FaceLighting = 'gouraud';
    s.AmbientStrength = 1;
    s.DiffuseStrength = 0.8;
    s.BackFaceLighting = 'lit';
    
    s.SpecularStrength = 1;
    s.SpecularColorReflectance = 1;
    s.SpecularExponent = 7;
    t = annotation('textbox','String','$F=6$, $H_-=0.28$ width is $0.18$ with periodic $y$--boundary condition','Interpreter','latex','EdgeColor','none');
    t.FontSize=20;
    t.Position=[0.43,0.85,0.4,0.12];
    axis off
    ax2=subplot(2,1,2);
    plot(shock(:,1),shock(:,2),'.')
    xlim([X_lower X_upper])
    ylim([0 0.18])
    t=xlabel('x','Interpreter','latex');
    t.FontSize=20;
    t=ylabel('y','Interpreter','latex');
    t.FontSize=15;
    formatt='t=%.1f';
    t=title(sprintf(formatt,(i-1)*dt),'Interpreter','latex');
    t.FontSize=20;
    mov(i)=getframe(f);
    clf(f)
end

v=VideoWriter(filenamev,'MPEG-4');
open(v)
writeVideo(v,mov)
close(v)