%% Example
% example section of the paper
% transaction process of the SEP mechanism
% version: matlab 2018a
clc;clear;
%% data
load=-3*[175.19,165.15,158.67,154.73,155.06,160.48,179.39,177.6,186.81,206.96,228.61,236.1,...
     242.18,243.6,248.86,255.79,256,246.74,245.97,237.35,237.31,232.67,195.93,195.6];
WindPower=3*[44,70.2,76,82,84,84,100,100,78,64,100,92,32,24,28,30,25,36,56,84,80,78,82,52];
PV=2*[0, 0, 0, 0, 0, 40, 76,82,95,76,85,100,108,112,106,95,85, 73, 32, 0,0,0,0,0];
load_p=load/sum(load);%shape of load curve

T=length(load);
rep=@(x)(repmat(x,1,T));
Pmax=[110;300;max(WindPower);max(PV)];Pmin=[55;200;0;0];%G1, G2, G3 and G4 respectively
Ramp=rep([100;300;max(WindPower);max(PV)]);
ae=0.02*[0.1;0.05;0.03;0.03];be=[65;60;46;50];ce=[0;0;0;0];%cost coefficients of substitute energy (G1-G4)
as=0.08*[0.02;0.01];bs=[35;30];cs=[0;0];%cost coefficients of regulation energy (G1-G2)
Sprice_pre=27.2;%preset price of regulation energy
Renew_Rate=(max(WindPower)+max(PV))/(sum(Pmax))
%% stage 1: substitute energy market
        Ng=size(Pmax,1);%G1-G4 participate in the substitute energy market
        Pe=sdpvar(Ng,T);%energy curve
        P_DS1=Pe(1,:)-load_p*sum(Pe(1,:));%regulation demand of energy curve
        P_DS2=Pe(2,:)-load_p*sum(Pe(2,:));
        P_DS3=Pe(3,:)-load_p*sum(Pe(3,:));
        P_DS4=Pe(4,:)-load_p*sum(Pe(4,:));
        
        f_TotalEnergyBalan=(sum(sum(Pe,1))+sum(load)==0);
        f_Gbound1=(Pe<=rep(Pmax))+(Pe>=rep(Pmin))+(Pe(3,:)<=WindPower)+(Pe(4,:)<=PV);
        f_Gramp1=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                 (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));
        f_stage1=f_TotalEnergyBalan+f_Gbound1+f_Gramp1;
        
        obj_energy=sum(ae.*(sum(Pe,2)).^2+be.*sum(Pe,2)+ce);%cost of substitute energy
        obj_RRM=Sprice_pre*(norm(P_DS1,1)+norm(P_DS2,1)+norm(P_DS3,1)+norm(P_DS4,1));%potential cost of regulation demand
        obj_stage1=obj_energy+obj_RRM;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage1,obj_stage1,os);
        
        Pe_opt=value(Pe);
        RR1=load_p*sum(Pe_opt(1,:));%substitute energy curve considering regulation responsibility
        RR2=load_p*sum(Pe_opt(2,:));
        RR3=load_p*sum(Pe_opt(3,:));
        RR4=load_p*sum(Pe_opt(4,:));
        AbandWP=WindPower-Pe_opt(3,:);
        AbandPV=PV-Pe_opt(4,:);
        
        P_DS_system=load+sum(Pe_opt,1);%system regulation demand
        SDopt=[];%regulation demand of each entity
        SDopt(1,:)=norm(Pe_opt(1,:)-load_p*sum(Pe_opt(1,:)),1);
        SDopt(2,:)=norm(Pe_opt(2,:)-load_p*sum(Pe_opt(2,:)),1);
        SDopt(3,:)=norm(Pe_opt(3,:)-load_p*sum(Pe_opt(3,:)),1);
        SDopt(4,:)=norm(Pe_opt(4,:)-load_p*sum(Pe_opt(4,:)),1);
        
        SEP_each=2*ae.*sum(Pe_opt,2)+be;
        Price_SubstituteEnergy=max(SEP_each)
        revenue_energy=Price_SubstituteEnergy.*sum(Pe_opt,2);%revenue of substitute energy of each entity
        
        % indirect payment
        Ng_ip=size(Pmax,1);
        Pe_ip=sdpvar(Ng_ip,T);
        P_DS1_ip=Pe_ip(1,:)-load_p*sum(Pe_ip(1,:));
        P_DS2_ip=Pe_ip(2,:)-load_p*sum(Pe_ip(2,:));
        P_DS3_ip=Pe_ip(3,:)-load_p*sum(Pe_ip(3,:));
        P_DS4_ip=Pe_ip(4,:)-load_p*sum(Pe_ip(4,:));
        f_TotalEnergyBalan_ip=(sum(sum(Pe_ip,1))+sum(load)==0);
        f_Gbound1_ip=(Pe_ip<=rep(Pmax))+(Pe_ip>=rep(Pmin))+(Pe_ip(3,:)<=WindPower)+(Pe_ip(4,:)<=PV);
        f_Gramp1_ip=(Pe_ip(:,2:end)-Pe_ip(:,1:end-1)>=-Ramp(:,2:end))+...
                    (Pe_ip(:,2:end)-Pe_ip(:,1:end-1)<= Ramp(:,2:end));
        f_stage1_ip=f_TotalEnergyBalan_ip+f_Gbound1_ip+f_Gramp1_ip;
        obj_energy_ip=sum(ae.*(sum(Pe_ip,2)).^2+be.*sum(Pe_ip,2)+ce);
        % Set potential cost of regulation demand as zero, to calculate indirect payment of each entity
        obj_RRM_ip=0*Sprice_pre*(norm(P_DS1_ip,1)+norm(P_DS2_ip,1)+norm(P_DS3_ip,1)+norm(P_DS4_ip,1));
        obj_stage1_ip=obj_energy_ip+obj_RRM_ip;
        os_ip=sdpsettings('solver','gurobi','verbose',0);
        sol_ip=solvesdp(f_stage1_ip,obj_stage1_ip,os_ip);
        Pe_opt_ip=value(Pe_ip);
        lose_SE_market_share=sum(Pe_opt,2)-sum(Pe_opt_ip,2);
        indirect_payment=Price_SubstituteEnergy.*lose_SE_market_share

%% stage 2: regulation energy market
as_ES=[0.0024];bs_ES=[40];cs_ES=[0];%cost coefficients of regulation energy (GES)
Rs_ESmax=2000;Rs_ESmin=400;
Ps_ESmax=200; Ps_ESmin=-200;
        Ps_G=sdpvar(2,T);%regulation curve (G1-G2)
        Ps_ES=sdpvar(1,T);%regulation curve (GES)
        Rs_ES=sdpvar(1,T);
        f_TimeEnergyBalan=(sum(Ps_G,1)+sum(Ps_ES,1)+P_DS_system==0);
        f_Gbound2=(Pe_opt(1:2,:)+Ps_G<=rep(Pmax(1:2)))+(Pe_opt(1:2,:)+Ps_G>=rep(Pmin(1:2)))+(sum(Ps_G,2)==0);
        f_Gramp2=(Pe_opt(1,2:end)+Ps_G(1,2:end)-Pe_opt(1,1:end-1)-Ps_G(1,1:end-1)>=-Ramp(1,2:end))+...
                 (Pe_opt(1,2:end)+Ps_G(1,2:end)-Pe_opt(1,1:end-1)-Ps_G(1,1:end-1)<= Ramp(1,2:end));
        f_ESbound=(Ps_ES<=Ps_ESmax)+(Ps_ES>=Ps_ESmin)+(Rs_ES(1,2:end)==Rs_ES(1,1:end-1)-1*Ps_ES(1,1:end-1))...
                 +(Rs_ES<=rep(Rs_ESmax))+(Rs_ES>=rep(Rs_ESmin))+(Rs_ES(1,1)==Rs_ESmax/2)+(sum(Ps_ES,2)==0);
        f_stage2=f_TimeEnergyBalan+f_Gbound2+f_Gramp2+f_ESbound;
        
        obj_Gregulation=sum(as.*sum(abs(Ps_G),2).^2  +bs.*sum(abs(Ps_G),2)  +cs);
        obj_ESregulation=sum(as_ES.*sum(abs(Ps_ES),2).^2  +bs_ES.*sum(abs(Ps_ES),2)  +cs_ES);
        obj_stage2=obj_Gregulation+obj_ESregulation;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage2,obj_stage2,os);
        
        Ps_Gopt=value(Ps_G);
        Ps_ESopt=value(Ps_ES);
        Rs_ESopt=value(Rs_ES);
        P_Gfinal(1:2,:)=Pe_opt(1:2,:)+Ps_Gopt;
        
        priceS_G=max(2*as.*sum(abs(Ps_Gopt),2)+bs);
        priceS_ES=2*as_ES.*sum(abs(Ps_ESopt),2)+bs_ES;
        Price_RegulationEnergy=max(priceS_G,priceS_ES)

        revenue_regulation_G1G2=[Price_RegulationEnergy*norm(Ps_Gopt(1,:),1);Price_RegulationEnergy*norm(Ps_Gopt(2,:),1)];
        revenue_regulation_ES=Price_RegulationEnergy*norm(Ps_ESopt,1);
        revenue_regulation_total=sum(revenue_regulation_G1G2)+revenue_regulation_ES;
        
        payment_load=-sum(revenue_energy)-revenue_regulation_total %payment of load
        RegulationEnergy_TransRate=revenue_regulation_total/sum(revenue_energy)
        
%% display
% Figure3 Energy curves supplied by market entities
% and corresponding substitute energy curves considering regulation responsibility
    dc=lines;
    G0=zeros(1,24);
    figure(3)
    x=1:1:T;
    y=G0;
    subplot(2,1,1);
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    % G1 energy curve
    y=Pe_opt(1,:)+G0;
    figure(3)
    subplot(2,1,1);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:))
    % G2 energy curve
    YY=[];
    XX=x;
    y=Pe_opt(1,:)+Pe_opt(2,:)+G0;
    YY=Pe_opt(1,:)+G0;
    figure(3)
    subplot(2,1,1);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(20,:))
    % G3 energy curve
    YY=[];
    XX=x;
    y=Pe_opt(1,:)+Pe_opt(2,:)+Pe_opt(3,:)+G0;
    YY=Pe_opt(1,:)+Pe_opt(2,:)+G0;
    figure(3)
    subplot(2,1,1);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(30,:));hold on;
    % G4 energy curve
    YY=[];
    XX=x;
    y=Pe_opt(1,:)+Pe_opt(2,:)+Pe_opt(3,:)+Pe_opt(4,:)+G0;
    YY=Pe_opt(1,:)+Pe_opt(2,:)+Pe_opt(3,:)+G0;
    figure(3)
    subplot(2,1,1);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(40,:));hold on;
    % G1 substitute energy curve considering regulation responsibility
    figure(3)
    x=1:1:T;
    y=G0;
    subplot(2,1,2);
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    y=RR1+G0;
    figure(3)
    subplot(2,1,2);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:))
    % G2 substitute energy curve considering regulation responsibility
    YY=[];
    XX=x;
    y=RR2+RR1+G0;
    YY=RR1+G0;
    figure(3)
    subplot(2,1,2);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(20,:))
    % G3 substitute energy curve considering regulation responsibility
    YY=[];
    XX=x;
    y=RR1+RR2+RR3+G0;
    YY=RR1+RR2+G0;
    figure(3)
    subplot(2,1,2);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(30,:));hold on;
    % G4 substitute energy curve considering regulation responsibility
    YY=[];
    XX=x;
    y=RR1+RR2+RR3+RR4+G0;
    YY=RR1+RR2+RR3+G0;
    figure(3)
    subplot(2,1,2);
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(40,:));hold on;

% Figure 4 System regulation demand and its components
    figure(4)
    x=1:1:T;
    y=G0;
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    y=P_DS_system;
    figure(4)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:));
    plot(Pe_opt(2,:)-load_p*sum(Pe_opt(2,:)));hold on;
    plot(Pe_opt(3,:)-load_p*sum(Pe_opt(3,:)));hold on;
    plot(Pe_opt(4,:)-load_p*sum(Pe_opt(4,:)));hold on;
    
% Figure 5 Regulation curves and contributions of regulation energy of market entities
    figure(5)
    x=1:1:T;
    y=RR1;
    subplot(5,1,1)
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    % G1
    y=Pe_opt(1,:);
    figure(5)
    subplot(3,1,1)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:))
    YY=[];
    XX=x;
    y=Ps_Gopt(1,:)+Pe_opt(1,:);
    YY=Pe_opt(1,:);
    figure(5)
    subplot(3,1,1)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(20,:))
    % G2
    figure(5)
    x=1:1:T;
    y=RR2;
    subplot(3,1,2)
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    y=Pe_opt(2,:);
    figure(5)
    subplot(3,1,2)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:))
    YY=[];
    XX=x;
    y=Ps_Gopt(2,:)+Pe_opt(2,:);
    YY=Pe_opt(2,:);
    figure(5)
    subplot(3,1,2)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(20,:))

    %GES
    figure(5)
    x=1:1:T;
    y=G0;
    subplot(3,1,3)
    plot(x,y);
    XX=x;
    YY=y;
    hold on;
    y=Ps_ESopt(1,:);
    figure(5)
    subplot(5,1,5)
    plot(x,y)
    hold on;
    XX=[XX x(end:-1:1)];
    YY=[YY y(end:-1:1)];
    patch(XX,YY,dc(10,:))
    
% Table 1 Revenues and payments of entities in SEP-based market mechanism
Table1=[revenue_energy(1),0,revenue_regulation_G1G2(1),0,revenue_energy(1)+revenue_regulation_G1G2(1);
        revenue_energy(2),0,revenue_regulation_G1G2(2),0,revenue_energy(2)+revenue_regulation_G1G2(2);
        revenue_energy(3),0,0,0,revenue_energy(3);
        revenue_energy(4),0,0,0,revenue_energy(4);
        0,0,revenue_regulation_ES,0,revenue_regulation_ES]