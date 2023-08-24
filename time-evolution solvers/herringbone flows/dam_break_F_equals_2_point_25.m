clc
clear
filepath=strcat(pwd,'\dam_break_F_equals_2_point_25');
N=5000;
mov(1:N)= struct('cdata',[],'colormap',[]);
filenamev='dambreak_F_equals_2_point_25';
X_lower=0.005;
X_upper=19.995;
dx=0.01;
Y_lower=0.005;
Y_upper=4.995;
dy=0.01;
[X,Y] = meshgrid(X_lower:dx:X_upper,Y_lower:dy:Y_upper);
[m,n]=size(X);
dt=0.2;
for i=1:N
    if i==1
        f=figure;
        f.Color = 'black';
    end
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

    ax=axes;
    s=surf(X,Y,H);
    s.EdgeColor='none';
    view(3)
    ax.XLim = [0 20];
    ax.ZLim = [0 5];
    ax.CameraPosition =[12 -6 5];
    ax.CameraTarget = [10 2.5 1];
    ax.CameraUpVector = [0 0 1];
    ax.CameraViewAngle = 25;
    ax.DataAspectRatio = [1 1 1];
    l1 = light;
    theta=60/180*pi;
    alpha=1/2*pi;
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
    t = annotation('textbox','String','$F=2.25$, $H_R=0.7$, channel width is $5$ and length is $20$','Interpreter','latex','EdgeColor','none');
    t.FontSize=20;
    t.Color='white';
    t.Position=[0.27,0.05,0.4,0.12];
    t.FitBoxToText='on';
    formatt='t=%.1f';
    t=annotation('textbox','String',sprintf(formatt,(i-1)*dt),'Interpreter','latex','FontSize',20,'EdgeColor','none');
    t.FontSize=20;
    t.Color='white';
    t.Position=[0.5,0.8,0.4,0.12];
    axis off
    mov(i)=getframe(f);
    clf(f)
end
v=VideoWriter(filenamev,'MPEG-4');
open(v)
writeVideo(v,mov)
close(v)
