%���ĸ���
%˫�㲩�ģ�KKT����
%κ�|, �«h, ����, et al. �������Ӳ��ĵ�����С�������̶��۲��Լ��綯����������[J]. ��������, 2015(4).
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
Ce=sdpvar(24,1);%���
z=binvar(24,1);%���۵�״̬
u=binvar(24,1);%����״̬
Pb=sdpvar(24,1);%��ǰ����
Pb_day=sdpvar(24,1);%ʵʱ����
Ps_day=sdpvar(24,1);%ʵʱ�۵�
Pdis=sdpvar(24,1);%���ܷŵ�
Pch=sdpvar(24,1);%���ܳ��
Pc1=sdpvar(24,1);%һ�೵��繦��
Pc2=sdpvar(24,1);%���೵��繦��
Pc3=sdpvar(24,1);%���೵��繦��
S=sdpvar(24,1);%��������
for t=2:24
    S(t)=S(t-1)+0.9*Pch(t)-Pdis(t)/0.9;
end
C=[lb<=Ce<=ub,mean(Ce)==0.5,Pb>=0,Ps_day<=Pdis,Pb_day>=0,Pb_day<=1000*z,Ps_day>=0,Ps_day<=1000*(1-z),Pch>=0,Pch<=1000*u,Pdis>=0,Pdis<=1000*(1-u)];%�߽�Լ��
C=[C,Pc1+Pc2+Pc3+Pch-Pdis==Pb+Pb_day-Ps_day];%����ƽ��
C=[C,sum(0.9*Pch-Pdis/0.9)==0,S(24)==2500,S>=0,S<=5000];%SOCԼ��
L_u=sdpvar(1,3);%���������ʽԼ�����������պ���
L_lb=sdpvar(24,3);%��繦������Լ�����������պ���
L_ub=sdpvar(24,3);%��繦������Լ�����������պ���
L_T=sdpvar(24,3);%������ʱ��Լ�����������պ���
f=50*L_u(1)*(0.9*24-9.6)+20*L_u(2)*(0.9*24-9.6)+10*L_u(3)*(0.9*24-9.6)+sum(sum(L_ub).*[50*3,20*3,10*3])+sum(price_s.*Ps_day-price_day_ahead.*Pb-price_b.*Pb_day);%Ŀ�꺯��
C=[C,Ce-L_u(1)*ones(24,1)-L_lb(:,1)-L_ub(:,1)-L_T(:,1)==0,Ce-L_u(2)*ones(24,1)-L_lb(:,2)-L_ub(:,2)-L_T(:,2)==0,Ce-L_u(3)*ones(24,1)-L_lb(:,3)-L_ub(:,3)-L_T(:,3)==0];%KKT����
C=[C,sum(Pc1)==50*(0.9*24-9.6),sum(Pc2)==20*(0.9*24-9.6),sum(Pc3)==10*(0.9*24-9.6)];%��������Լ��
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
b_lb=binvar(24,3);%��繦������Լ�����ɳڱ���
b_ub=binvar(24,3);%��繦������Լ�����ɳڱ���
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
