function [A, B, C] = GetMPCControlMatrix(VehConf, V)
    Cr = VehConf.cr;
    Cf = VehConf.cf;
    m  = VehConf.m;
    Lr = VehConf.b;
    Lf = VehConf.a;
    Iz = VehConf.Iz;
    A = zeros(6,6);
    A(1,2) = 1;
    
    A(2,2) = -(Cf+Cr)/(m*V);
    A(2,3) = (Cf+Cr)/m;
    A(2,4) = (-Cf*Lf+Cr*Lr)/(m*V);
    
    A(3,4) = 1;
    
    A(4,2) = -(Cf*Lf - Cr*Lr)/(Iz*V);
    A(4,3) =  (Cf*Lf - Cr*Lr)/Iz;
    A(4,4) = -(Cf*Lf^2 + Cr*Lr^2)/(Iz*V);
    
    A(5,6) = 1;
    
    B = zeros(6,2);
    B(2,1) = Cf/m; 
    B(4,1) =Lf*Cf/Iz;
    B(6,2) = -1;
    
    C = zeros(6,1);
    C(2) = (Cr*Lr-Cf*Lf)/m/V - V;
    C(4) =-(Cr*Lr^2+Cf*Lf^2)/Iz/V;
    C(6) = 1;
end
