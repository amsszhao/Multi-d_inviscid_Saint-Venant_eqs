clc
clear
load('sub_analytic.mat');
load('parameterspace.mat');
j=1;
filename=strcat('job',string(j+1),'.mat');
b=cell(500,1);
poolobj=parpool('local',64);
for i=1:1:500
    i
    b{i}=extract_sub(parameter(i+j*500,1),parameter(i+j*500,2),10^-12,10^-12,10^-14,firstcolumn_func,eigen_mat,eigen_vector,[0:0.2:6],[0.1],[30],[200],[300],[1000]);
    save(filename,'b');
end
delete(poolobj)
