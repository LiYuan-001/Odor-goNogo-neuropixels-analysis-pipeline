function [Nlx_TS,Nlx_O,P_task] = cal_TS_TP_onebox(Nlx_t_type,signal_ss)
%%% odor A
Nlx_TS.Hit1  = Nlx_t_type(signal_ss(:,2)==1 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);
Nlx_O.Hit1   =  signal_ss(signal_ss(:,2)==1 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);

Nlx_TS.Miss1 = Nlx_t_type(signal_ss(:,2)==1 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);
Nlx_O.Miss1  =  signal_ss(signal_ss(:,2)==1 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);

Nlx_TS.Im1   = Nlx_t_type(signal_ss(:,2)==1 & isnan(signal_ss(:,4)),1);
Nlx_O.Im1    =  signal_ss(signal_ss(:,2)==1 & isnan(signal_ss(:,4)),1);

%%% odor B
Nlx_TS.CR2   = Nlx_t_type(signal_ss(:,2)==2 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);
Nlx_O.CR2    =  signal_ss(signal_ss(:,2)==2 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);

Nlx_TS.FA2   = Nlx_t_type(signal_ss(:,2)==2 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);
Nlx_O.FA2    =  signal_ss(signal_ss(:,2)==2 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);

Nlx_TS.Im1   = Nlx_t_type(signal_ss(:,2)==2 & isnan(signal_ss(:,4)),1);
Nlx_O.Im1    =  signal_ss(signal_ss(:,2)==2 & isnan(signal_ss(:,4)),1);

%%% odor C
Nlx_TS.Hit3  = Nlx_t_type(signal_ss(:,2)==3 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);
Nlx_O.Hit3   =  signal_ss(signal_ss(:,2)==3 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);

Nlx_TS.Miss3 = Nlx_t_type(signal_ss(:,2)==3 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);
Nlx_O.Miss3  =  signal_ss(signal_ss(:,2)==3 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);

Nlx_TS.Im3   = Nlx_t_type(signal_ss(:,2)==3 & isnan(signal_ss(:,4)),1);
Nlx_O.Im3    =  signal_ss(signal_ss(:,2)==3 & isnan(signal_ss(:,4)),1);

%%% odor D
Nlx_TS.CR4  = Nlx_t_type(signal_ss(:,2)==4 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);
Nlx_O.CR4   =  signal_ss(signal_ss(:,2)==4 & signal_ss(:,4)==1 & signal_ss(:,5) == 0,1);

Nlx_TS.FA4  = Nlx_t_type(signal_ss(:,2)==4 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);
Nlx_O.FA4   =  signal_ss(signal_ss(:,2)==4 & signal_ss(:,4)==0 & signal_ss(:,5) == 0,1);

Nlx_TS.Im4   = Nlx_t_type(signal_ss(:,2)==4 & isnan(signal_ss(:,4)),1);
Nlx_O.Im4    =  signal_ss(signal_ss(:,2)==4 & isnan(signal_ss(:,4)),1);


%% Calculate the performance of task
P_task.A=sum(signal_ss(:,2)==1 & signal_ss(:,4)==1 & signal_ss(:,5) == 0)/sum(signal_ss(:,2)==1 & ~isnan(signal_ss(:,4)) & signal_ss(:,5) == 0)*100;
P_task.B=sum(signal_ss(:,2)==2 & signal_ss(:,4)==1 & signal_ss(:,5) == 0)/sum(signal_ss(:,2)==2 & ~isnan(signal_ss(:,4)) & signal_ss(:,5) == 0)*100;
P_task.C=sum(signal_ss(:,2)==3 & signal_ss(:,4)==1 & signal_ss(:,5) == 0)/sum(signal_ss(:,2)==3 & ~isnan(signal_ss(:,4)) & signal_ss(:,5) == 0)*100;
P_task.D=sum(signal_ss(:,2)==4 & signal_ss(:,4)==1 & signal_ss(:,5) == 0)/sum(signal_ss(:,2)==4 & ~isnan(signal_ss(:,4)) & signal_ss(:,5) == 0)*100;
end

