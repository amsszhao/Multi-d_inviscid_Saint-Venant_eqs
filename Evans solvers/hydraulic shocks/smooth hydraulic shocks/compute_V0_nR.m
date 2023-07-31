%compute the expansion of the nonradial solutions at +infty for eta=0
f=obj.F;
hr=obj.HR;

eigen=-(f*(16*hr^2*lam - 3*f^2*hr + 8*hr^(3/2)*lam + 8*hr^(5/2)*lam + 4*f^2 + f^2*hr^2 - 4*f^2*hr^(1/2) + 2*f^2*hr^(3/2) + 4*hr^2*lam^2 + 4*hr^3*lam^2 + 8*hr^(5/2)*lam^2 - 4*f^2*hr^(3/2)*lam - 4*f^2*hr*lam)^(1/2) - 2*f^2*hr + 2*f^2 + f^2*hr^(1/2) - f^2*hr^(3/2) + f*hr^(1/2)*(16*hr^2*lam - 3*f^2*hr + 8*hr^(3/2)*lam + 8*hr^(5/2)*lam + 4*f^2 + f^2*hr^2 - 4*f^2*hr^(1/2) + 2*f^2*hr^(3/2) + 4*hr^2*lam^2 + 4*hr^3*lam^2 + 8*hr^(5/2)*lam^2 - 4*f^2*hr^(3/2)*lam - 4*f^2*hr*lam)^(1/2) + 2*f^2*hr^(1/2)*lam + 2*f^2*hr*lam)/(2*(2*hr^2 + hr^(3/2) + hr^(5/2) - f^2*hr^(1/2)));
V=zeros(3,n);
V(:,1)=[(2*f^2*hr^(1/2)*lam*(hr^(1/2) + 1))/(f*(16*hr^2*lam - 3*f^2*hr + 8*hr^(3/2)*lam + 8*hr^(5/2)*lam + 4*f^2 + f^2*hr^2 - 4*f^2*hr^(1/2) + 2*f^2*hr^(3/2) + 4*hr^2*lam^2 + 4*hr^3*lam^2 + 8*hr^(5/2)*lam^2 - 4*f^2*hr^(3/2)*lam - 4*f^2*hr*lam)^(1/2) - f^2*hr + 2*f^2 - f^2*hr^(1/2) + 2*f^2*hr^(3/2)*lam + 2*f^2*hr*lam);1;0];
V(:,1)=V(:,1)/norm(V(:,1));
for i=2:n
    RHS=zeros(3,1);
    for j=3:min(6,i)
        RHS=RHS-(i+1-j)*obj.LR(j)*V(:,i+2-j);
    end
    for j=2:min(7,i)
        RHS=RHS+lam*obj.RR2{j}*V(:,i+1-j);
    end
    for j=2:min(6,i)
        RHS=RHS+obj.RR1{j}*V(:,i+1-j)+eigen*obj.RR4(j)*V(:,i+1-j);
    end
    for j3=1:min(i,7)
        for j1=1:min(i-1,i+1-j3)
            for j2=max(1,4-j1-j3):min(i-1,i+2-j1-j3)
                j4=i+3-j1-j2-j3;
                mat=lam*obj.RR2{j3};
                if j3<=6
                    mat=mat+obj.RR1{j3}+eigen*obj.RR4(j3)*eye(3);
                end
                RHS=RHS-V(:,j1)*V(:,j2)'*mat*V(:,j4);
            end
        end
    end
    V(:,i)=((i-1)*obj.LR(2)*eye(3)-(eye(3)-V(:,1)*V(:,1)')*(obj.RR1{1}+lam*obj.RR2{1}+eigen*obj.RR4(1)*eye(3)))^-1*RHS;
end


