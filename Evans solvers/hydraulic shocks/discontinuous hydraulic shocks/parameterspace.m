clc
clear
hr=0.01:0.01:0.99;
n=length(hr);
k=1;
for i=1:1:n
    v=(1/hr(i))^0.5;
    F_min=1/v+1/v^2;
    F_max=(hr(i)^0.5+1)/(hr(i)*(hr(i)^0.5-hr(i)+1))^0.5;
    f=0.02;
    while f<F_min
        f=f+0.02;
    end
    while f<F_max
        parameter(k,:)=[f hr(i)];
        k=k+1;
        f=f+0.02;
    end
end
