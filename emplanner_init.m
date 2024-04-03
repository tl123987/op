
%%%参数设置%%%%%%
DEG2RAD=pi/180;
RAD2DEG=180/pi;
%%%%整车参数%%%%%
cf=-155016;
cr=-92204;
m=1422;
Iz=1536.0;
la=1.225;
lb=2.68-1.225;
%%%%%%%横向LQR参数
LQR_Q1=50;
LQR_Q2=15;
LQR_Q3=30;
LQR_Q4=15;
LQR_R=10;

k=zeros(5000,4);
vx_break_point=zeros(1,5000);
for i=1:5000
    vx_break_point(i)=0.01*i;
    
    A=[0,1,0,0;
        0,(cf+cr)/(m*vx_break_point(i)),-(cf+cr)/m,(la*cf-lb*cr)/(m*vx_break_point(i));
        0,0,0,1;
        0,(la*cf-lb*cr)/(Iz*vx_break_point(i)),-(la*cf-lb*cr)/Iz,(la*la*cf+lb*lb*cr)/(Iz*vx_break_point(i))];
    B=[0;
        -cf/m;
        0;
        -la*cf/Iz];
LQR_Q=1*[LQR_Q1,0,0,0;
        0,LQR_Q2,0,0;
        0,0,LQR_Q3,0;
        0,0,0,LQR_Q4];
   k(i,:)=lqr(A,B,LQR_Q,LQR_R);
end
LQR_K1=k(:,1)';
LQR_K2=k(:,2)';
LQR_K3=k(:,3)';
LQR_K4=k(:,4)';
%%%%车辆初始位置
host_x_init=0; 
host_y_init=0;


