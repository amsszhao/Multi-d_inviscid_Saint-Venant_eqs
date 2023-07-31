%compute the expansion of the solutions at -infty for eta=0
L2=obj.L(2);
L3=obj.L(3);
L4=obj.L(4);
L5=obj.L(5);
L6=obj.L(6);

R11=obj.R1{1};
R12=obj.R1{2};
R13=obj.R1{3};
R14=obj.R1{4};
R15=obj.R1{5};
R16=obj.R1{6};

R21=obj.R2{1};
R22=obj.R2{2};
R23=obj.R2{3};
R24=obj.R2{4};
R25=obj.R2{5};
R26=obj.R2{6};
R27=obj.R2{7};

R41=obj.R4(1);
R42=obj.R4(2);
R43=obj.R4(3);
R44=obj.R4(4);
R45=obj.R4(5);
R46=obj.R4(6);

f=obj.F;
hr=obj.HR;
eigen=(f*(hr^(1/2)+1)*((8*lam+8*hr*lam-3*f^2*hr+4*hr*lam^2+16*hr^(1/2)*lam+f^2+4*f^2*hr^2+2*f^2*hr^(1/2)-4*f^2*hr^(3/2)+4*lam^2+8*hr^(1/2)*lam^2-4*f^2*hr^(3/2)*lam-4*f^2*hr*lam)^(1/2)-f*hr^(1/2)-f+2*f*hr+2*f*hr*lam))/(2*hr+4*hr^(1/2)-2*f^2*hr^2+2);
V1=[f-2*f*lam+f*hr^(1/2)+(8*lam+8*hr*lam-3*f^2*hr+4*hr*lam^2+16*hr^(1/2)*lam+f^2+4*f^2*hr^2+2*f^2*hr^(1/2)-4*f^2*hr^(3/2)+4*lam^2+8*hr^(1/2)*lam^2-4*f^2*hr^(3/2)*lam-4*f^2*hr*lam)^(1/2)-2*f*hr-2*f*hr^(1/2)*lam;2*f*lam*(hr^(1/2)+1);0];

V2=(1*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(+(lam*R22+R12+eigen*R42*eye(3))*V1);
V3=(2*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-L3*V2+(lam*R22+R12+eigen*R42*eye(3))*V2+(lam*R23+R13+eigen*R43*eye(3))*V1);
V4=(3*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-2*L3*V3-L4*V2+(lam*R22+R12+eigen*R42*eye(3))*V3+(lam*R23+R13+eigen*R43*eye(3))*V2+(lam*R24+R14+eigen*R44*eye(3))*V1);
V5=(4*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-3*L3*V4-2*L4*V3-L5*V2+(lam*R22+R12+eigen*R42*eye(3))*V4+(lam*R23+R13+eigen*R43*eye(3))*V3+(lam*R24+R14+eigen*R44*eye(3))*V2+(lam*R25+R15+eigen*R45*eye(3))*V1);
V6=(5*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-4*L3*V5-3*L4*V4-2*L5*V3-L6*V2+(lam*R22+R12+eigen*R42*eye(3))*V5+(lam*R23+R13+eigen*R43*eye(3))*V4+(lam*R24+R14+eigen*R44*eye(3))*V3+(lam*R25+R15+eigen*R45*eye(3))*V2+(lam*R26+R16+eigen*R46*eye(3))*V1);
V7=(6*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-5*L3*V6-4*L4*V5-3*L5*V4-2*L6*V3+(lam*R22+R12+eigen*R42*eye(3))*V6+(lam*R23+R13+eigen*R43*eye(3))*V5+(lam*R24+R14+eigen*R44*eye(3))*V4+(lam*R25+R15+eigen*R45*eye(3))*V3+(lam*R26+R16+eigen*R46*eye(3))*V2+(lam*R27)*V1);
V8=(7*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-6*L3*V7-5*L4*V6-4*L5*V5-3*L6*V4+(lam*R22+R12+eigen*R42*eye(3))*V7+(lam*R23+R13+eigen*R43*eye(3))*V6+(lam*R24+R14+eigen*R44*eye(3))*V5+(lam*R25+R15+eigen*R45*eye(3))*V4+(lam*R26+R16+eigen*R46*eye(3))*V3+(lam*R27)*V2);
V9=(8*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-7*L3*V8-6*L4*V7-5*L5*V6-4*L6*V5+(lam*R22+R12+eigen*R42*eye(3))*V8+(lam*R23+R13+eigen*R43*eye(3))*V7+(lam*R24+R14+eigen*R44*eye(3))*V6+(lam*R25+R15+eigen*R45*eye(3))*V5+(lam*R26+R16+eigen*R46*eye(3))*V4+(lam*R27)*V3);
V10=(9*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-8*L3*V9-7*L4*V8-6*L5*V7-5*L6*V6+(lam*R22+R12+eigen*R42*eye(3))*V9+(lam*R23+R13+eigen*R43*eye(3))*V8+(lam*R24+R14+eigen*R44*eye(3))*V7+(lam*R25+R15+eigen*R45*eye(3))*V6+(lam*R26+R16+eigen*R46*eye(3))*V5+(lam*R27)*V4);
V11=(10*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-9*L3*V10-8*L4*V9-7*L5*V8-6*L6*V7+(lam*R22+R12+eigen*R42*eye(3))*V10+(lam*R23+R13+eigen*R43*eye(3))*V9+(lam*R24+R14+eigen*R44*eye(3))*V8+(lam*R25+R15+eigen*R45*eye(3))*V7+(lam*R26+R16+eigen*R46*eye(3))*V6+(lam*R27)*V5);
V12=(11*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-10*L3*V11-9*L4*V10-8*L5*V9-7*L6*V8+(lam*R22+R12+eigen*R42*eye(3))*V11+(lam*R23+R13+eigen*R43*eye(3))*V10+(lam*R24+R14+eigen*R44*eye(3))*V9+(lam*R25+R15+eigen*R45*eye(3))*V8+(lam*R26+R16+eigen*R46*eye(3))*V7+(lam*R27)*V6);
V13=(12*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-11*L3*V12-10*L4*V11-9*L5*V10-8*L6*V9+(lam*R22+R12+eigen*R42*eye(3))*V12+(lam*R23+R13+eigen*R43*eye(3))*V11+(lam*R24+R14+eigen*R44*eye(3))*V10+(lam*R25+R15+eigen*R45*eye(3))*V9+(lam*R26+R16+eigen*R46*eye(3))*V8+(lam*R27)*V7);
V14=(13*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-12*L3*V13-11*L4*V12-10*L5*V11-9*L6*V10+(lam*R22+R12+eigen*R42*eye(3))*V13+(lam*R23+R13+eigen*R43*eye(3))*V12+(lam*R24+R14+eigen*R44*eye(3))*V11+(lam*R25+R15+eigen*R45*eye(3))*V10+(lam*R26+R16+eigen*R46*eye(3))*V9+(lam*R27)*V8);
V15=(14*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-13*L3*V14-12*L4*V13-11*L5*V12-10*L6*V11+(lam*R22+R12+eigen*R42*eye(3))*V14+(lam*R23+R13+eigen*R43*eye(3))*V13+(lam*R24+R14+eigen*R44*eye(3))*V12+(lam*R25+R15+eigen*R45*eye(3))*V11+(lam*R26+R16+eigen*R46*eye(3))*V10+(lam*R27)*V9);
V16=(15*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-14*L3*V15-13*L4*V14-12*L5*V13-11*L6*V12+(lam*R22+R12+eigen*R42*eye(3))*V15+(lam*R23+R13+eigen*R43*eye(3))*V14+(lam*R24+R14+eigen*R44*eye(3))*V13+(lam*R25+R15+eigen*R45*eye(3))*V12+(lam*R26+R16+eigen*R46*eye(3))*V11+(lam*R27)*V10);
V17=(16*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-15*L3*V16-14*L4*V15-13*L5*V14-12*L6*V13+(lam*R22+R12+eigen*R42*eye(3))*V16+(lam*R23+R13+eigen*R43*eye(3))*V15+(lam*R24+R14+eigen*R44*eye(3))*V14+(lam*R25+R15+eigen*R45*eye(3))*V13+(lam*R26+R16+eigen*R46*eye(3))*V12+(lam*R27)*V11);
V18=(17*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-16*L3*V17-15*L4*V16-14*L5*V15-13*L6*V14+(lam*R22+R12+eigen*R42*eye(3))*V17+(lam*R23+R13+eigen*R43*eye(3))*V16+(lam*R24+R14+eigen*R44*eye(3))*V15+(lam*R25+R15+eigen*R45*eye(3))*V14+(lam*R26+R16+eigen*R46*eye(3))*V13+(lam*R27)*V12);
V19=(18*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-17*L3*V18-16*L4*V17-15*L5*V16-14*L6*V15+(lam*R22+R12+eigen*R42*eye(3))*V18+(lam*R23+R13+eigen*R43*eye(3))*V17+(lam*R24+R14+eigen*R44*eye(3))*V16+(lam*R25+R15+eigen*R45*eye(3))*V15+(lam*R26+R16+eigen*R46*eye(3))*V14+(lam*R27)*V13);
V20=(19*L2*eye(3)-(R11+lam*R21+eigen*R41*eye(3)))^-1*(-18*L3*V19-17*L4*V18-16*L5*V17-15*L6*V16+(lam*R22+R12+eigen*R42*eye(3))*V19+(lam*R23+R13+eigen*R43*eye(3))*V18+(lam*R24+R14+eigen*R44*eye(3))*V17+(lam*R25+R15+eigen*R45*eye(3))*V16+(lam*R26+R16+eigen*R46*eye(3))*V15+(lam*R27)*V14);

V=[V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12 V13 V14 V15 V16 V17 V18 V19 V20];