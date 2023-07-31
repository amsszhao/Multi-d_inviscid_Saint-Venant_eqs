%compute the expansion of the nonradial solutions at +infty for eta~=0
f=obj.F;
hr=obj.HR;

eigen_matR_v=obj.eigen_matR(f,hr,eta,lam);
eigen_vectorR_v=obj.eigen_vectorR(f,hr,eta,lam);
[~,index]=min(real([eigen_matR_v(1,1) eigen_matR_v(2,2) eigen_matR_v(3,3)]));
eigen=eigen_matR_v(index,index);
V=zeros(3,n);
V(:,1)=eigen_vectorR_v(:,index);
V(:,1)=V(:,1)/norm(V(:,1));

for i=2:n
    RHS=zeros(3,1);
    for j=3:min(6,i)
        RHS=RHS-(i+1-j)*obj.LR(j)*V(:,i+2-j);
    end
    for j=2:min(7,i)
        RHS=RHS+lam*obj.RR2{j}*V(:,i+1-j)+eta*obj.RR3{j}*V(:,i+1-j);
    end
    for j=2:min(6,i)
        RHS=RHS+obj.RR1{j}*V(:,i+1-j)+eigen*obj.RR4(j)*V(:,i+1-j);
    end
    for j3=1:min(7,i)
        for j1=1:min(i-1,i+1-j3)
            for j2=max(1,4-j1-j3):min(i-1,i+2-j1-j3)
                j4=i+3-j1-j2-j3;
                mat=lam*obj.RR2{j3}+eta*obj.RR3{j3};
                if j3<=6
                    mat=mat+obj.RR1{j3}+eigen*obj.RR4(j3)*eye(3);
                end
                RHS=RHS-V(:,j1)*V(:,j2)'*mat*V(:,j4);
            end
        end
    end
    V(:,i)=((i-1)*obj.LR(2)*eye(3)-(eye(3)-V(:,1)*V(:,1)')*(obj.RR1{1}+lam*obj.RR2{1}+eta*obj.RR3{1}+eigen*obj.RR4(1)*eye(3)))^-1*RHS;
end


