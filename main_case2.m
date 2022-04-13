%% case 2
% case2 reflects the cost causality of regulation demand.
% assume that there are only two wind power generators (denoted as WPG1 and WPG2) and one GES in a 100% RE power system
% all of which are price makers.
% version: matlab 2018a
%% case 2-1 (LMP mechanism)
        clc;clear;
        load=-3*[175.19,165.15,158.67,154.73,155.06,160.48,179.39,177.6,186.81,206.96,228.61,236.1,...
             242.18,243.6,248.86,255.79,256,246.74,245.97,237.35,237.31,232.67,195.93,195.6];
        WindPower1=-0.55*load;
        WindPower2=-0.55*load;
        load_p=load/sum(load);
        T=length(load);
        rep=@(x)(repmat(x,1,T));
        Pmax=[max(WindPower1);max(WindPower2)];Pmin=[0;0];
        Ramp=rep([max(WindPower1);max(WindPower2)]);
        a=[0.03;0.03];b=[46;46];c=[0;0];

        Ng=size(Pmin,1);
        Pe=sdpvar(Ng,T);
        as_ES=[0];bs_ES=[40];cs_ES=[0];
        Rs_ESmax=10000;Rs_ESmin=2000;
        Ps_ESmax=1000; Ps_ESmin=-1000;
        Ps_ES=sdpvar(1,T);
        Rs_ES=sdpvar(1,T);
        
        f_Balan=(sum(Pe,1)+Ps_ES+load==0);
        f_Gbound=(Pe>=rep(Pmin))+(Pe(1,:)<=WindPower1)+(Pe(2,:)<=WindPower2);
        f_Gramp=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));    
        f_ESbound=(Ps_ES<=Ps_ESmax)+(Ps_ES>=Ps_ESmin)+(Rs_ES(1,2:end)==Rs_ES(1,1:end-1)-1*Ps_ES(1,1:end-1))...
                 +(Rs_ES<=rep(Rs_ESmax))+(Rs_ES>=rep(Rs_ESmin))+(Rs_ES(1,1)==Rs_ESmax/2)+(sum(Ps_ES,2)==0);
        f_all=f_Balan+f_Gbound+f_Gramp+f_ESbound;
        
        obj_G=sum(sum(rep(a).*Pe.^2+rep(b).*Pe+rep(c)));
        obj_ES=sum(sum(rep(as_ES(1)).*abs(Ps_ES).^2+rep(bs_ES(1)).*abs(Ps_ES)+rep(cs_ES(1))));
        obj_all=obj_G+obj_ES;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_all,obj_all,os);
        
        Pe_opt=value(Pe);
        Ps_ESopt=value(Ps_ES);
        AbandWP1=WindPower1-Pe_opt(1,:);
        AbandWP2=WindPower2-Pe_opt(2,:);
        LMP=-dual(f_Balan)
        revenueWP=sum(LMP.*Pe_opt,2);
        revenueES=sum(LMP.*Ps_ESopt);
        payment_load=-sum(revenueWP)-revenueES;
        result_LMP=[revenueWP;revenueES;payment_load];
        
%% case 2-2 (LMP mechanism)
        WindPower=3*[44,70.2,76,82,84,84,100,100,78,64,100,92,32,24,28,30,25,36,56,84,80,78,82,52];
        WindPower1=-0.55*load;
        WindPower2=-0.55*WindPower/sum(WindPower)*sum(load);
        Pmax=[max(WindPower1);max(WindPower2)];Pmin=[0;0];
        Ramp=rep([max(WindPower1);max(WindPower2)]);
        Ng=size(Pmin,1);
        Pe=sdpvar(Ng,T);
        Ps_ES=sdpvar(1,T);
        Rs_ES=sdpvar(1,T);
        f_Balan=(sum(Pe,1)+Ps_ES+load==0);
        f_Gbound=(Pe>=rep(Pmin))+(Pe(1,:)<=WindPower1)+(Pe(2,:)<=WindPower2);
        f_Gramp=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                 (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));
        f_ESbound=(Ps_ES<=Ps_ESmax)+(Ps_ES>=Ps_ESmin)+(Rs_ES(1,2:end)==Rs_ES(1,1:end-1)-1*Ps_ES(1,1:end-1))...
                 +(Rs_ES<=rep(Rs_ESmax))+(Rs_ES>=rep(Rs_ESmin))+(Rs_ES(1,1)==Rs_ESmax/2)+(sum(Ps_ES,2)==0);
        f_all=f_Balan+f_Gbound+f_Gramp+f_ESbound;
        
        obj_G=sum(sum(rep(a).*Pe.^2+rep(b).*Pe+rep(c)));
        obj_ES=sum(sum(rep(as_ES(1)).*abs(Ps_ES).^2+rep(bs_ES(1)).*abs(Ps_ES)+rep(cs_ES(1))));
        obj_all=obj_G+obj_ES;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_all,obj_all,os);
        
        Pe_opt=value(Pe);
        Ps_ESopt=value(Ps_ES);
        AbandWP1=WindPower1-Pe_opt(1,:);
        AbandWP2=WindPower2-Pe_opt(2,:);
        LMP=-dual(f_Balan)
        revenueWP=sum(LMP.*Pe_opt,2);
        revenueES=sum(LMP.*Ps_ESopt);
        payment_load=-sum(revenueWP)-revenueES;
        
        result_LMP(1:2,2)=revenueWP;
        result_LMP(3,2)=revenueES;
        result_LMP(4,2)=payment_load;
   
%% case 2-1 (SEP mechanism)
        WindPower1=-0.55*load;
        WindPower2=-0.55*load;
        Pmax=[max(WindPower1);max(WindPower2)];Pmin=[0;0];
        Ramp=rep([max(WindPower1);max(WindPower2)]);
        ae=0.0013;%Seting ae=0.0013 is to ensure that load payments under the two mechanisms are equal in case2-1 for comparison and contrast.
        be=b;ce=c;
        Sprice_pre=bs_ES;
        
        %substitute energy market
        Ng=size(Pmax,1);
        Pe=sdpvar(Ng,T);
        P_DS1=Pe(1,:)-load_p*sum(Pe(1,:));
        P_DS2=Pe(2,:)-load_p*sum(Pe(2,:));
 
        f_TotalEnergyBalan=(sum(sum(Pe,1))+sum(load)==0);
        f_Gbound1=(Pe<=rep(Pmax))+(Pe>=rep(Pmin))+(Pe(1,:)<=WindPower1)+(Pe(2,:)<=WindPower2);
        f_Gramp1=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                 (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));
        f_stage1=f_TotalEnergyBalan+f_Gbound1+f_Gramp1;
        
        obj_energy=sum(ae.*(sum(Pe,2)).^2+be.*sum(Pe,2)+ce);
        obj_RRM=Sprice_pre*(norm(P_DS1,1)+norm(P_DS2,1));
        obj_stage1=obj_energy+obj_RRM;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage1,obj_stage1,os);
        
        Pe_opt=value(Pe);
        AbandWP1=WindPower1-Pe_opt(1,:);
        AbandWP2=WindPower2-Pe_opt(2,:);
        SDopt=[];
        SDopt(1,:)=norm(Pe_opt(1,:)-load_p*sum(Pe_opt(1,:)),1);
        SDopt(2,:)=norm(Pe_opt(2,:)-load_p*sum(Pe_opt(2,:)),1);
        P_DS_system=load+sum(Pe_opt,1);
        
        SEP_each=2*ae.*sum(Pe_opt,2)+be;
        Price_SubstituteEnergy=max(SEP_each)
        revenueWP=Price_SubstituteEnergy.*sum(Pe_opt,2);
       
        %regulation energy market
        Ps_ES=sdpvar(1,T);
        Rs_ES=sdpvar(1,T);
        f_TimeEnergyBalan=(Ps_ES+P_DS_system==0);
        f_ESbound=(Ps_ES<=Ps_ESmax)+(Ps_ES>=Ps_ESmin)+(Rs_ES(1,2:end)==Rs_ES(1,1:end-1)-1*Ps_ES(1,1:end-1))...
                 +(Rs_ES<=rep(Rs_ESmax))+(Rs_ES>=rep(Rs_ESmin))+(Rs_ES(1,1)==Rs_ESmax/2)+(sum(Ps_ES,2)==0);
        f_stage2=f_TimeEnergyBalan+f_ESbound;
        
        obj_ESregulation=sum(as_ES.*sum(abs(Ps_ES),2).^2  +bs_ES.*sum(abs(Ps_ES),2)  +cs_ES);
        obj_stage2=obj_ESregulation;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage2,obj_stage2,os);
                
        Ps_ESopt=value(Ps_ES);
        Price_RegulationEnergy=2*as_ES.*sum(abs(Ps_ESopt),2)+bs_ES;
        revenueES=Price_RegulationEnergy*norm(Ps_ESopt,1);
        payment_load=-sum(revenueWP)-revenueES;
        result_SEP=[revenueWP;revenueES;payment_load];

%% case 2-2 (SEP mechanism)
        WindPower=3*[44,70.2,76,82,84,84,100,100,78,64,100,92,32,24,28,30,25,36,56,84,80,78,82,52];
        WindPower1=-0.55*load;
        WindPower2=-0.55*WindPower/sum(WindPower)*sum(load);
        Pmax=[max(WindPower1);max(WindPower2)];Pmin=[0;0];
        Ramp=rep([max(WindPower1);max(WindPower2)]);
        %substitute energy market
        Ng=size(Pmax,1);
        Pe=sdpvar(Ng,T);
        P_DS1=Pe(1,:)-load_p*sum(Pe(1,:));
        P_DS2=Pe(2,:)-load_p*sum(Pe(2,:));
        f_TotalEnergyBalan=(sum(sum(Pe,1))+sum(load)==0);
        f_Gbound1=(Pe<=rep(Pmax))+(Pe>=rep(Pmin))+(Pe(1,:)<=WindPower1)+(Pe(2,:)<=WindPower2);
        f_Gramp1=(Pe(:,2:end)-Pe(:,1:end-1)>=-Ramp(:,2:end))+...
                 (Pe(:,2:end)-Pe(:,1:end-1)<= Ramp(:,2:end));
        f_stage1=f_TotalEnergyBalan+f_Gbound1+f_Gramp1;
        obj_energy=sum(ae.*(sum(Pe,2)).^2+be.*sum(Pe,2)+ce);
        obj_RRM=Sprice_pre*(norm(P_DS1,1)+norm(P_DS2,1));
        obj_stage1=obj_energy+obj_RRM;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage1,obj_stage1,os);
        Pe_opt=value(Pe);
        AbandWP1=WindPower1-Pe_opt(1,:);
        AbandWP2=WindPower2-Pe_opt(2,:);
        SDopt=[];
        SDopt(1,:)=norm(Pe_opt(1,:)-load_p*sum(Pe_opt(1,:)),1);
        SDopt(2,:)=norm(Pe_opt(2,:)-load_p*sum(Pe_opt(2,:)),1);    
        P_DS_system=load+sum(Pe_opt,1);
        
        SEP_each=2*ae.*sum(Pe_opt,2)+be;
        Price_SubstituteEnergy=max(SEP_each)
        revenueWP=Price_SubstituteEnergy.*sum(Pe_opt,2);
        
        %regulation energy market
        Ps_ES=sdpvar(1,T);
        Rs_ES=sdpvar(1,T);
        f_TimeEnergyBalan=(Ps_ES+P_DS_system==0);
        f_ESbound=(Ps_ES<=Ps_ESmax)+(Ps_ES>=Ps_ESmin)+(Rs_ES(1,2:end)==Rs_ES(1,1:end-1)-1*Ps_ES(1,1:end-1))...
                 +(Rs_ES<=rep(Rs_ESmax))+(Rs_ES>=rep(Rs_ESmin))+(Rs_ES(1,1)==Rs_ESmax/2)+(sum(Ps_ES,2)==0);
        f_stage2=f_TimeEnergyBalan+f_ESbound;
        
        obj_ESregulation=sum(as_ES.*sum(abs(Ps_ES),2).^2  +bs_ES.*sum(abs(Ps_ES),2)  +cs_ES);
        obj_stage2=obj_ESregulation;
        os=sdpsettings('solver','gurobi','verbose',0);
        sol=solvesdp(f_stage2,obj_stage2,os);
        
        Ps_ESopt=value(Ps_ES);
        Price_RegulationEnergy=2*as_ES.*sum(abs(Ps_ESopt),2)+bs_ES;
        revenueES=Price_RegulationEnergy*norm(Ps_ESopt,1);
        revenueWP=revenueWP;
        payment_load=-sum(revenueWP)-revenueES;
        
        result_SEP(1:2,2)=revenueWP;
        result_SEP(3,2)=revenueES;
        result_SEP(4,2)=payment_load;

%% display
% Figure 7 Comparison of SEP mechanism and LMP mechanism in case 2-1 and case 2-2
        figure(7)
        subplot(2,2,1);
        bar(result_LMP/(10e+3));
        subplot(2,2,2);
        bar(result_SEP/(10e+3));
        subplot(2,2,3);
        plot(-load);hold on;
        plot(WindPower1);hold on;%WPG1
        plot(WindPower1);hold on;%Case2-1:WPG2
        plot(WindPower2);        %Case2-2:WPG2
        result_LMP
        result_SEP
        