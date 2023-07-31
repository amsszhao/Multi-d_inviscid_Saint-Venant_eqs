%The following function computes expansions of the solutions near the sonic
%point
function W_con=compute_W(n,lam,eta,F,L,R1,R2,R3)
W_con=zeros(3,2,n+1);
W_con(:,:,1)=[(F*lam*(F-2))/(F+1)^2,(eta*3i)/(F^2*(F-2));0,(4*F-3*lam+3*F*lam+2*F^2+2)/(F^2*(F-2));((F-2)*(4*F-3*lam+3*F*lam+2*F^2+2))/(6*(F+1)),0];
W_con(1:2,:,2)=(R1(1:2,:,1)+lam*R2(1:2,:,1)+eta*R3(1:2,:,1))*W_con(:,:,1);
for i=2:n
    RHS=zeros(3,2);
    for j=2:min(4,i)
        RHS=RHS+R1(:,:,j)*W_con(:,:,i+1-j);
    end
    for j=2:min(5,i)
        RHS=RHS+lam*R2(:,:,j)*W_con(:,:,i+1-j);
    end
    for j=2:min(6,i)
        RHS=RHS+eta*R3(:,:,j)*W_con(:,:,i+1-j);
    end
    for j=3:min(6,i)
        RHS=RHS-(i+1-j)*L(:,:,j)*W_con(:,:,i+2-j);
    end
    W_con(3,:,i)=((i-1)*L(3,3,2)-R1(3,3,1)-R2(3,3,1)*lam-R3(3,3,1)*eta)^-1*(-(i-1)*L(3,1:2,2)*W_con(1:2,:,i)+(R1(3,1:2,1)+R2(3,1:2,1)*lam+R3(3,1:2,1)*eta)*W_con(1:2,:,i)+RHS(3,:));
    RHS=RHS-(i-1)*L(:,:,2)*W_con(:,:,i)+(R1(:,:,1)+R2(:,:,1)*lam+R3(:,:,1)*eta)*W_con(:,:,i);
    W_con(1:2,:,i+1)=1/i*RHS(1:2,:);
end
end
