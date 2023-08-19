clc
clear

filepath1=strcat(pwd,'\dam_break_F_equals_2_point_13');
filepath2=strcat(pwd,'\dam_break_F_equals_2_point_14');
filepath3=strcat(pwd,'\dam_break_F_equals_2_point_15');

N=7000;
mov(1:N)= struct('cdata',[],'colormap',[]);
filenamev='dam_break_comparision';
X_lower=0.005;
X_upper=19.995;
dx=0.01;
Y_lower=0.005;
Y_upper=4.995;

dy=0.01;
[X,Y] = meshgrid(X_lower:dx:X_upper,Y_lower:dy:Y_upper);
[m,n]=size(X);

dt=0.1;

for i=1:N
    if i==1
       f=figure;
       f.Color='black';
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
    ax=subplot(3,1,1);
    filename1=fullfile(filepath1,sprintf(format,i-1));
    fileID=fopen(filename1);
    C=textscan(fileID,'%f %f %f','HeaderLines',9);
    fclose(fileID);
    H=C{1};
    H=reshape(H,[n,m]);
    H=transpose(H);
    s=surf(X,Y,H);
    s.EdgeColor='none';
    view(3)
    ax.XLim = [0 20];
    ax.ZLim = [0 5];
    ax.CameraPosition =[10.5 -6 5];
    ax.CameraTarget = [10.5 2.5 1];
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
     t = annotation('textbox','String','$F=2.13$, $H_R=0.7$, channel width is $5$ and length is $20$','Interpreter','latex','EdgeColor','none');
    t.FontSize=16;
    t.Color='white';
    t.Position=[0.31,0.6,0.4,0.12];
    formatt='t=%.1f';
    t=annotation('textbox','String',sprintf(formatt,(i-1)*dt),'Interpreter','latex','FontSize',20,'EdgeColor','none');
    t.FontSize=25;
    t.Color='white';
    t.Position=[0.46,0.85,0.4,0.12];
    t.FitBoxToText='on';
    axis off
    ax=subplot(3,1,2);
    filename2=fullfile(filepath2,sprintf(format,i-1));
    fileID=fopen(filename2);
    C=textscan(fileID,'%f %f %f','HeaderLines',9);
    fclose(fileID);
    H=C{1};
    H=reshape(H,[n,m]);
    H=transpose(H);
    s=surf(X,Y,H);
    s.EdgeColor='none';
    view(3)
    ax.XLim = [0 20];
    ax.ZLim = [0 5];
    ax.CameraPosition =[10.5 -6 5];
    ax.CameraTarget = [10.5 2.5 1];
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
     t = annotation('textbox','String','$F=2.14$, $H_R=0.7$, channel width is $5$ and length is $20$','Interpreter','latex','EdgeColor','none');
    t.FontSize=16;
    t.Color='white';
    t.Position=[0.31,0.3,0.4,0.12];
    t.FitBoxToText='on';
    axis off
    ax=subplot(3,1,3);
    filename3=fullfile(filepath3,sprintf(format,i-1));
    fileID=fopen(filename3);
    C=textscan(fileID,'%f %f %f','HeaderLines',9);
    fclose(fileID);
    H=C{1};
    H=reshape(H,[n,m]);
    H=transpose(H);
    s=surf(X,Y,H);
    s.EdgeColor='none';
    view(3)
    ax.XLim = [0 20];
    ax.ZLim = [0 5];
    ax.CameraPosition =[10.5 -6 5];
    ax.CameraTarget = [10.5 2.5 1];
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
     t = annotation('textbox','String','$F=2.15$, $H_R=0.7$, channel width is $5$ and length is $20$','Interpreter','latex','EdgeColor','none');
    t.FontSize=16;
    t.Color='white';
    t.FitBoxToText='on';
    t.Position=[0.31,0,0.4,0.12];
    axis off
    mov(i)=getframe(f);
    clf(f)
end
v=VideoWriter(filenamev,'MPEG-4');
open(v)
writeVideo(v,mov)
close(v)
