clc
clear
hr=0.01:0.01:0.99;
n=length(hr);
k=1;
for i=1:1:n
    v=(1/hr(i))^0.5;
    F_max=1/v+1/v^2;
    f=0.02;
    while f<F_max
        parameter(k,:)=[f hr(i)];
        f=f+0.02;
        k=k+1;
    end
end
