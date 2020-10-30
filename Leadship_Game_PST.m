%论文复现
%双层博弈，KKT条件
%魏韡, 陈玥, 刘锋, et al. 基于主从博弈的智能小区代理商定价策略及电动汽车充电管理[J]. 电网技术, 2015(4).
clear
clc
price_day_ahead=[0.35;0.33;0.3;0.33;0.36;0.4;0.44;0.46;0.52;0.58;0.66;0.75;0.81;0.76;0.8;0.83;0.81;0.75;0.64;0.55;0.53;0.47;0.40;0.37];
price_b=1.2*price_day_ahead;
price_s=1.2*price_day_ahead;
lb=0.8*price_day_ahead;
ub=1.2*price_day_ahead;
T_1=[1;1;1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;1;1;1];
T_2=[1;1;1;1;1;1;1;1;0;0;0;0;1;1;1;0;0;0;0;1;1;1;1;1];
T_3=[0;0;0;0;0;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;0;0;0];
Ce=sdpvar(24,1);%电价
z=binvar(24,1);%购售电状态
u=binvar(24,1);%储能状态
Pb=sdpvar(24,1);%日前购电
Pb_day=sdpvar(24,1);%实时购电
Ps_day=sdpvar(24,1);%实时售电
Pdis=sdpvar(24,1);%储能放电
Pch=sdpvar(24,1);%储能充电
Pc1=sdpvar(24,1);%一类车充电功率
Pc2=sdpvar(24,1);%二类车充电功率
Pc3=sdpvar(24,1);%三类车充电功率
S=sdpvar(24,1);%储荷容量
for t=2:24
    S(t)=S(t-1)+0.9*Pch(t)-Pdis(t)/0.9;
end
C=[lb<=Ce<=ub,mean(Ce)==0.5,Pb>=0,Ps_day<=Pdis,Pb_day>=0,Pb_day<=1000*z,Ps_day>=0,Ps_day<=1000*(1-z),Pch>=0,Pch<=1000*u,Pdis>=0,Pdis<=1000*(1-u)];%边界约束
C=[C,Pc1+Pc2+Pc3+Pch-Pdis==Pb+Pb_day-Ps_day];%能量平衡
C=[C,sum(0.9*Pch-Pdis/0.9)==0,S(24)==2500,S>=0,S<=5000];%SOC约束
L_u=sdpvar(1,3);%电量需求等式约束的拉格朗日函数
L_lb=sdpvar(24,3);%充电功率下限约束的拉格朗日函数
L_ub=sdpvar(24,3);%充电功率上限约束的拉格朗日函数
L_T=sdpvar(24,3);%充电可用时间约束的拉格朗日函数
f=50*L_u(1)*(0.9*24-9.6)+20*L_u(2)*(0.9*24-9.6)+10*L_u(3)*(0.9*24-9.6)+sum(sum(L_ub).*[50*3,20*3,10*3])+sum(price_s.*Ps_day-price_day_ahead.*Pb-price_b.*Pb_day);%目标函数
C=[C,Ce-L_u(1)*ones(24,1)-L_lb(:,1)-L_ub(:,1)-L_T(:,1)==0,Ce-L_u(2)*ones(24,1)-L_lb(:,2)-L_ub(:,2)-L_T(:,2)==0,Ce-L_u(3)*ones(24,1)-L_lb(:,3)-L_ub(:,3)-L_T(:,3)==0];%KKT条件
C=[C,sum(Pc1)==50*(0.9*24-9.6),sum(Pc2)==20*(0.9*24-9.6),sum(Pc3)==10*(0.9*24-9.6)];%电量需求约束
for t=1:24
    if T_1(t)==0
        C=[C,Pc1(t)==0];
    else
        C=[C,L_T(t,1)==0];
    end
    if T_2(t)==0
        C=[C,Pc2(t)==0];
    else
        C=[C,L_T(t,2)==0];
    end
    if T_3(t)==0
        C=[C,Pc3(t)==0];
    else
        C=[C,L_T(t,3)==0];
    end
end
b_lb=binvar(24,3);%充电功率下限约束的松弛变量
b_ub=binvar(24,3);%充电功率上限约束的松弛变量
M=1000;
for t=1:24
    if T_1(t)==0
        C=[C,L_ub(t,1)==0,b_ub(t,1)==1,b_lb(t,1)==1];
    else
        C=[C,L_lb(t,1)>=0,L_lb(t,1)<=M*b_lb(t,1),Pc1(t)>=0,Pc1(t)<=M*(1-b_lb(t,1)),Pc1(t)<=50*3,50*3-Pc1(t)<=M*b_ub(t,1),L_ub(t,1)<=0,L_ub(t,1)>=M*(b_ub(t,1)-1)];
    end
    if T_2(t)==0
        C=[C,L_ub(t,2)==0,b_ub(t,2)==1,b_lb(t,2)==1];
    else
        C=[C,L_lb(t,2)>=0,L_lb(t,2)<=M*b_lb(t,2),Pc2(t)>=0,Pc2(t)<=M*(1-b_lb(t,2)),Pc2(t)<=20*3,20*3-Pc2(t)<=M*b_ub(t,2),L_ub(t,2)<=0,L_ub(t,2)>=M*(b_ub(t,2)-1)];
    end
    if T_3(t)==0
        C=[C,L_ub(t,3)==0,b_ub(t,3)==1,b_lb(t,3)==1];
    else
        C=[C,L_lb(t,3)>=0,L_lb(t,3)<=M*b_lb(t,3),Pc3(t)>=0,Pc3(t)<=M*(1-b_lb(t,3)),Pc3(t)<=10*3,10*3-Pc3(t)<=M*b_ub(t,3),L_ub(t,3)<=0,L_ub(t,3)>=M*(b_ub(t,3)-1)];
    end
end
ops=sdpsettings('solver','cplex');
solvesdp(C,-f,ops);
Pc=[double(Pc1),double(Pc2),double(Pc3)];
Pb=double(Pb);
Ps_day=double(Ps_day);
Pb_day=double(Pb_day);
S=double(S);
Pch=double(Pch);
Pdis=double(Pdis);
Cost_total=double(f)
Price_Charge=double(Ce);
clear ans b_lb b_ub C Ce f L_lb L_ub L_T L_u lb M ops Pc1 Pc2 Pc3 price_b price_day_ahead price_s t T_1 T_2 T_3 u ub z
