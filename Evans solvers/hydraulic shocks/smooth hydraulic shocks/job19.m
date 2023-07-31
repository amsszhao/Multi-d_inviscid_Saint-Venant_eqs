clc
clear
load('smooth_analytic.mat');
load('parameterspace.mat');
j=18;
filename=strcat('job',string(j+1),'.mat');
b=cell(250,1);
poolobj=parpool('local',64);
for i=1:1:250
    i
    b{i}=extract_smooth(parameter(i+j*250,1),parameter(i+j*250,2),10^-12,10^-12,10^-14,eigen_mat,eigen_matR,eigen_vector,eigen_vectorR,[0:0.2:6],[0.1],[30],[200],[200],[1000]);
    save(filename,'b');
end
delete(poolobj)