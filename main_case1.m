%% case 1
% case 1 tests the performance of valuing regulation contribution.
% version: matlab 2018a
clc;clear;
%% data
load0=-3*[175.19,165.15,158.67,154.73,155.06,160.48,179.39,177.6,186.81,206.96,228.61,236.1,...
     242.18,243.6,248.86,255.79,256,246.74,245.97,237.35,237.31,232.67,195.93,195.6];
load_p=load0/sum(load0);
WindPower=3*[44,70.2,76,82,84,84,100,100,78,64,100,92,32,24,28,30,25,36,56,84,80,78,82,52];
PV=2*[0, 0, 0, 0, 0, 40, 76,82,95,76,85,100,108,112,106,95,85, 73, 32, 0,0,0,0,0];
T=length(load0);
rep=@(x)(repmat(x,1,T));
%% LMP mechanism
% The two LMPs (in the 11th and 13th time period)
% will guide a GES to provide regulation contribution for the system
% and the GES will arbitrage on the price difference.
Pmax=[300;300;max(WindPower);max(PV)];Pmin=[150;200;0;0];
Ramp=rep([300;300;max(WindPower);max(PV)]);
a=[0.1;0.05;0.03;0.03];b=[65;60;46;50];c=[0;0;0;0];
Renew_Rate=max(WindPower+PV)/(sum(Pmax));
q=1;
ES_contrib=[];
revenueES_LMP=[];
ES_power=zeros(1,T);
for i=0:0.25:9
        ES_power(1,11)=-10*i;
        ES_power(1,13)= 10*i;
        load=load0+ES_power;
        Ng=size(Pmax,1);
        Pe=sdpvar(Ng,T);
        
        f_EnergyBalan=(sum(Pe,1)+load==0);
        f_Gbound=(Pe<=rep(Pmax))+(Pe>=rep(Pmin))+(Pe(3,:)<=WindPower)+(Pe(4,:)<=PV);
        f_Gramp=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));        
        f_all=f_EnergyBalan+f_Gbound+f_Gramp;
        
        obj_LMP=sum(sum(rep(a).*Pe.^2+rep(b).*Pe+rep(c)));
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_all,obj_LMP,os);
        
        LMP=-dual(f_EnergyBalan);
        LMP11(q)=LMP(11);
        LMP13(q)=LMP(13);
        revenueES_LMP(q)=sum(LMP.*ES_power);
        ES_contrib(q)=norm(ES_power,1)/2;
        q=q+1;
end
%Figure 6 Comparison of SEP mechanism and LMP mechanism in case 1
figure(6)
subplot(2,1,1);
yyaxis left;
bar(ES_contrib,revenueES_LMP);
yyaxis right;
plot(ES_contrib,LMP11);hold on;
plot(ES_contrib,LMP13)

%% SEP mechanism
load=-3*[175.19,165.15,158.67,154.73,155.06,160.48,179.39,177.6,186.81,206.96,228.61,236.1,...
     242.18,243.6,248.86,255.79,256,246.74,245.97,237.35,237.31,232.67,195.93,195.6];
load_p=load/sum(load);
q=1;
ae=0.02*a;be=b;ce=c;
Sprice_pre=6.9;
revenueES_SEP=[];
Price_RegulationEnergy=[];

    %substitute energy market
        Ng=size(Pmax,1);
        Pe=sdpvar(Ng,T);
        P_DS1=Pe(1,:)-load_p*sum(Pe(1,:));
        P_DS2=Pe(2,:)-load_p*sum(Pe(2,:));
        P_DS3=Pe(3,:)-load_p*sum(Pe(3,:));
        P_DS4=Pe(4,:)-load_p*sum(Pe(4,:));
        
        f_TotalEnergyBalan=(sum(sum(Pe,1))+sum(load)==0);
        f_Gbound1=(Pe<=rep(Pmax))+(Pe>=rep(Pmin))+(Pe(3,:)<=WindPower)+(Pe(4,:)<=PV);
        f_Gramp1=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                 (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));
        f_stage1=f_TotalEnergyBalan+f_Gbound1+f_Gramp1;
        
        obj_energy=sum(ae.*(sum(Pe,2)).^2+be.*sum(Pe,2)+ce);
        obj_RRM=Sprice_pre*(norm(P_DS1,1)+norm(P_DS2,1)+norm(P_DS3,1)+norm(P_DS4,1));
        obj_stage1=obj_energy+obj_RRM;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage1,obj_stage1,os);
        
        Pe_opt=value(Pe);
        RR1=load_p*sum(Pe_opt(1,:));
        RR2=load_p*sum(Pe_opt(2,:));
        RR3=load_p*sum(Pe_opt(3,:));
        RR4=load_p*sum(Pe_opt(4,:));
        P_DS_system=load+sum(Pe_opt,1);
        
  for i=0:0.25:9    
    %regulation energy market
        ES_power=zeros(1,T);
        ES_power(1,11)=-10*i;
        ES_power(1,13)= 10*i;
        as_ES=[0.01;0];bs_ES=[9;10];cs_ES=[0;0];
        Ps_ES=sdpvar(2,T);
        f_TimeEnergyBalan=(sum(Ps_ES,1)+P_DS_system==0);
        f_ESbound=(sum(Ps_ES,2)==[0;0])+(Ps_ES(1,:)==ES_power);%simplified
        f_stage2=f_TimeEnergyBalan+f_ESbound;
        obj_ESregulation=sum(as_ES.*sum(abs(Ps_ES),2).^2  +bs_ES.*sum(abs(Ps_ES),2)  +cs_ES);
        obj_stage2=obj_ESregulation;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage2,obj_stage2,os);

        Ps_ESopt=value(Ps_ES);
        priceS_ES=2*as_ES.*sum(abs(Ps_ESopt),2)+bs_ES;
        Price_RegulationEnergy(q)=max(priceS_ES);
        revenueES_SEP(q)=Price_RegulationEnergy(q)*norm(Ps_ESopt(1,:),1);
        ES_contrib(q)=norm(ES_power,1)/2;
        q=q+1;
        
end
figure(6)
subplot(2,1,2);
yyaxis left;
bar(ES_contrib,revenueES_SEP);
yyaxis right;
plot(ES_contrib,Price_RegulationEnergy);

