function [sys,x0,str,ts] = lateral_control_mpc2(t,x,u,flag)
% 该程序功能：横纵向控制（适配gen）
% 程序版本 V1.0，MATLAB版本：R2022a,采用S函数的标准形式
% 程序编写日期 2022.10.26
% 最近一次改写 2022.10.27
% u 对应mux的顺序 u(1)=X0 u(2)=Y0 误差改版
% 状态量=[ed,ed_dot,ephi,ephi_dot]，状态量是横向误差，横向误差导数，航向角误差，航向角误差导数
% 控制量是前轮转角
switch flag
 case 0
  [sys,x0,str,ts] = mdlInitializeSizes; % Initialization
  
 case 2
  sys = mdlUpdates(t,x,u); 
  
 case 3
  sys = mdlOutputs(t,x,u); 
 

 case {1,4,9} 
  sys = [];
  
 otherwise
  error(['unhandled flag = ',num2str(flag)]); 
end

function [sys,x0,str,ts] = mdlInitializeSizes

sizes = simsizes;
sizes.NumContStates  = 0;
sizes.NumDiscStates  = 4;%这里应该是状态空间变量
sizes.NumOutputs     = 2;%包括控制量和其他检测量
sizes.NumInputs      = 15;%输入量个数
sizes.DirFeedthrough = 1; %直接贯通标识，不知道有什么用，默认值
sizes.NumSampleTimes = 1;%采样时间个数
sys = simsizes(sizes); %默认值

x0=[0;0;0;0];%状态量初始化
% x0=[0;0;0;0;0;0;0]
global U;
 U=0;
str = [];  %默认值           
ts  = [0.01 0];  %采样周期，偏移量  0.01 

function sys = mdlUpdates(t,x,u)
sys = x;
function sys = mdlOutputs(~,x,u)

%%
%%%%建立模型
global U ;
U=0;
%车辆参数
u(15);%后轮
u(14);%前轮
VehConf.cf=-148000;
% VehConf.cf=u(14);
VehConf.cr=-82204;
% VehConf.cr=u(15);
VehConf.Iz=1536.7;
VehConf.m=1412;
VehConf.a=1.22;
VehConf.b=2.680-1.22;
%%
tic
Nx=4;
Nu=1;
Np=20;
Nc=5;
Ny=4;

%%输入转换
vx=u(1)/3.6+0.1;
vy=u(2)/3.6;
yo=u(3);
xo=u(4);
xr=u(5);
yr=u(6);
theta=u(7);
yaw=u(8)*3.14/180;%弧度
wr=u(9)*3.14/180;%弧度
dbate=u(10)*3.14/180;
bate=u(11)*3.14/180;
kr=u(12);
vp=u(13)/3.6;
vact=sqrt(vx^2+vy^2);

kesi=zeros(Nx+Nu,1);
U=zeros(1,1);
%%误差向量
kesi(1)=(yo-yr)*cos(yaw-theta)-(xo-xr)*sin(yaw-theta);%横向误差
% kesi(2)=vact*sin(yaw-theta);
kesi(2)=vx*sin(yaw-theta)+vy*cos(yaw-theta);
kesi(3)=yaw-theta;%航向误差
kesi(4)=wr-vx*kr; %
kesi(5)=U(1,1);

dT=0.05;%采样时间0.01
[A1, B1, G1] = getcontrol(VehConf, vx);
A2 = eye(4) + dT*A1;
B2 = dT*B1;
G2=dT*G1;
A_cell=cell(2,2);
B_cell=cell(2,1);
A_cell{1,1}=A2;%6*6
A_cell{1,2}=B2;%6*1
A_cell{2,1}=zeros(Nu,Nx);
A_cell{2,2}=eye(Nu);
B_cell{1,1}=B2;
B_cell{2,1}=eye(Nu);
A=cell2mat(A_cell);%8X8
B=cell2mat(B_cell);%8X1
%整合
G_cell=cell(2,1);
G_cell{1,1}=G2;
G_cell{2,1}=zeros(Nu,1);
G=cell2mat(G_cell);%8X1

Q_cell=cell(Np,Np);
for i=1:1:Np
    for j=1:1:Np
        if i==j
             Q_cell{i,j}=[1000000 0 0 0;0 1000000 0 0;0 0 1000000 0;0 0 0 1000000];%10000000
        else 
           Q_cell{i,j}=zeros(Ny,Ny);               
        end
    end 
end
Q=cell2mat(Q_cell);%矩阵合并 80X80
R=2000000*eye(Nu*Nc);%nc=10,nu=1 1X5
C=[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0];%4*5
PSI_cell=cell(Np,1);%状态量前面的系数
%PSI为20*1,结果是A的i次方取前4行的4*5矩阵，预测时域为20，PSI为20个4*5的矩阵纵着放
for i=1:1:Np
    PSI_cell{i,1}=C*A^i;
end
PSI=cell2mat(PSI_cell);
%%%%%%%%%%%%%%%%20个A(k+1,2...)参数
thetar_cell=cell(Np,Nc);%20*10
for i=1:1:Np
    for j=1:1:Nc
        if j<=i
            thetar_cell{i,j}=C*A^(i-j)*B;
        else
            thetar_cell{i,j}=zeros(Ny,Nu);
        end
    end
end
thetar=cell2mat(thetar_cell);
tao_cell=cell(Np,Np);%20*20
PHI_cell=cell(Np,1);%20*1
for p=1:1:Np
    PHI_cell{p,1}=kr*(wr-vx*kr);%yaw_dot=kr*s_dot
    for q=1:1:Np
        if q<=p
        tao_cell{p,q}=C*A^(p-q)*G;%4*5 5*5 5*1
        else 
        tao_cell{p,q}=zeros(Ny,Nu);%4*1
        end 
    end
end
YITA_REF_CELL=[0;0;0;0];%质心侧偏角
%张量计算
Yita_ref=kron(ones(Np,1),YITA_REF_CELL);%%%%%%
PHI=cell2mat(PHI_cell);
tao=cell2mat(tao_cell);
error=PSI*kesi+tao*PHI;%应该是20个预测序列4*20列
H1=(thetar'*Q*thetar+R); %1X10
H=(H1'+H1)/2;
f=2*error'*Q*thetar;
%%  约束生成区域
%控制量约束
umin=-0.8;
umax=0.8;
Umin=kron(ones(Nc,1),umin);
Umax=kron(ones(Nc,1),umax);
A_t=zeros(Nc,Nc);
for i=1:1:Nc
    for j=1:1:Nc
        if j<=i
        A_t(i,j)=1;
        else
        A_t(i,j)=0;
        end
    end
end
A_I=kron(A_t,eye(Nu));
U_t=kron(ones(Nc,1),U);
delta_umin=-0.6;
delta_umax=0.6;

A_cons_cell=cell(2,1);
A_cons_cell{1,1}=A_I;
A_cons_cell{2,1}=-A_I;
A=cell2mat(A_cons_cell);
b_cons_cell=cell(2,1);
b_cons_cell{1,1}=Umax-U_t;
b_cons_cell{2,1}=-Umin+U_t;
b=cell2mat(b_cons_cell);

% ycmax=[-5;-5;-3;-5]; % 横摆角和纵向位置的约束
% ycmin=[5;5;3;5];
% Ycmax=kron(ones(Np,1),ycmax);
% Ycmin=kron(ones(Np,1),ycmin);
% A_cons={A_I ;-A_I;thetar;-thetar};
% b_cons={Umax-U_t;-Umin+U_t;Ycmax-PSI*kesi-tao*PHI;-Ycmin+PSI*kesi+tao*PHI};
% 
% A=cell2mat(A_cons);%（求解方程）状态量不等式约束增益矩阵，转换为绝对值的取值范围
% b=cell2mat(b_cons);%（求解方程）状态量不等式约束的取值


lb=delta_umin;
ub=delta_umax;
%%  开始求解过程，使用二次规划求解
options = optimset('Algorithm','active-set');
% options = optimset('Algorithm','interior-point');
x_start=eye(Nc,1)*dT;
[X,fval,exitflag]=quadprog(H,f,A,b,[],[],lb,ub,x_start,options);

sys = [X(1),kesi(1)];


