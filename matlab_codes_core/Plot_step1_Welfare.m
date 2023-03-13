%PLOTGAINS
%  
%   Current Policy との厚生比較 
%

% close all;
% clear all;
clc ;

%厚生の出力：全体、勤労、引退
agent = 'all';
% agent = 'worker'
% agent = 'retiree'

type = 'linear';
windowSize = 7; % 移動平均の期間

% scenarioを選びます
Scenario = 1;  % Baseline =1, retirement age ext=2, tfp high growth=3

%% policy 2から6のグラフを全部出します
for x = 2:4
 Policy   = x;  % Current =1, RR50=2, Med copay 10%up =3, LTC copay 10%up =4 
                % RR50 &  Med copay 10% down = 5, RR50 &  LTC copay 10% down = 6

%%
folder_name={'../file_out/grid1000_New/'} ;
 
  if Policy ==1       %Current
    if Scenario ==1
%          newpolicy = 'RR50';                range = [-40,20]; ee=-36;
    elseif   Scenario ==2
         newpolicy='Retire_ext'; range = [-5,2]; ee=1;
    elseif   Scenario ==3
         newpolicy= 'high_TFP'; range = [-5,10]; ee=9;
    end
 end   

 if Policy ==2       %RR50=2
    if Scenario ==1
         newpolicy = 'RR50';               range = [-5,4]; ee=3.5;
    elseif   Scenario ==2
         newpolicy='Retire_ext_RR50'; range = [-5,10]; ee=9;
    elseif   Scenario ==3
         newpolicy= 'high_TFP_RR50'; range = [-5,10]; ee=9;
    end
 end   
        
 if Policy ==3       %Med copay 10%up
    if Scenario ==1
         newpolicy = 'Med10';         range = [-5,4]; ee=3.5;
    elseif   Scenario ==2
         newpolicy='Retire_ext_MED10';  range = [-5,10]; ee=9;
    elseif   Scenario ==3
         newpolicy= 'high_TFP_MED10'; range = [-5,10]; ee=9;
    end
 end   
 
  if Policy ==4       % LTC copay 20%up
    if Scenario ==1
         newpolicy = 'Ltc10';        range = [-5,4]; ee=3.5;
    elseif   Scenario ==2
         newpolicy='Retire_ext_LTC10';  range = [-5,10]; ee=9; %range = [-5,4]; ee=3.5;
    elseif   Scenario ==3
         newpolicy= 'high_TFP_Ltc10'; range = [-5,10]; ee=9;
    end
  end   
  
  if Policy ==5       % RR50 & MED copay 10% Down
    if Scenario ==1
         newpolicy = 'RR50_MED10';   range = [-7.5,5]; ee=4;
    elseif   Scenario ==2
         newpolicy='Retire_ext_RR50_MED10';  range = [-5,10]; ee=9;
    elseif   Scenario ==3
         newpolicy= 'high_TFP_RR50_MED10'; range = [-5,10]; ee=9;
    end
 end   
 
  if Policy ==6       % RR50 & LTC copay 10% Down
    if Scenario ==1
         newpolicy = 'RR50_LTC10';   range = [-7.5,5]; ee=4;
    elseif   Scenario ==2
         newpolicy='Retire_ext_RR50_LTC10';  range = [-5,10]; ee=9;
    elseif   Scenario ==3
         newpolicy= 'high_TFP_RR50_Ltc10'; range = [-5,10]; ee=9;
    end
 end   
        
 %% new policy
 switch newpolicy
   %% Scenario 1 
 case 'Med10'
title_g= { '% Welfare change MED copay 30% & Cur. Policy in Baseline Scenario'};
  load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_MED_up.dat'] )
   newpolt  =trantf_muprop_sec2_tfp1_gz07_gn_grate3_MED_up;  
   
 case 'Med_Ltc30'
title_g= { '% Welfare changes Med and LTC copay 30% vs Current Policy'};
load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3.dat'] )
newpolt =trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3;

 case 'Ltc10'
title_g= { '% Welfare change LTC copay 30% & Cur. Policy in Baseline Scenario'};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_LTC_up.dat'] )
   newpolt  =trantf_muprop_sec2_tfp1_gz07_gn_grate3_LTC_up;  
  
  case 'RR50'
 title_g= { '%  Welfare change RR50% & Cur. Policy in Baseline Scenario '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3;
 
 case 'RR50_LTC10'
 title_g= { '%  Welfare change RR50% & LTC copay 10% down & Cur. Policy in Baseline Scenario '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_LTC.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_LTC;
 
 case 'RR50_MED10'
 title_g= { '%  Welfare change RR50% & MED copay 10% down & Cur. Policy in Baseline Scenario '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_MED.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_MED;
   
   
      %% Scenario 2
   case 'Retire_ext'; %
     title_g= { '% Welfare retire extension vs Current Policy'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire; 
     
  case 'Retire_ext_Med_Ltc30'
     title_g= { '% Welfare changes Med and LTC copay 30% vs Current Policy'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_retire.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_retire;   
    
  case 'Retire_ext_RR50'; %
     title_g= { '% Welfare changes RR 50%  & Cur. Policy in Retire Ext. Scenario '};   
      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire; 
     
 case 'Retire_ext_MED10'; %
     title_g= { '% Welfare MED Copay 30% & Cur. Policy in Retire Ext. Scenario'};   
      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_MED_up.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_MED_up;       
     
 case 'Retire_ext_LTC10'; %
     title_g= { '% Welfare LTC Copay 30% & Cur. Policy in Retire Ext. Scenario'};   
      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_LTC_up.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_LTC_up;     
     
%    case 'Retire_ext_RR50'
%  title_g= { '%  Welfare change RR50% & Cur. Policy in Retire Ext. Scenario '};
%  load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire.dat'] )
%  newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire;
 
 case 'Retire_ext_RR50_LTC10'
 title_g= { '%  Welfare change RR50% & LTC copay 10% down & Cur. Policy in Retire Ext. Scenario '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire_LTC.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire_LTC;
 
 case 'Retire_ext_RR50_MED10'
 title_g= { '%  Welfare change RR50% & MED copay 10% down & Cur. Policy in Retire Ext. Scenario '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire_MED.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire_MED;
  
     %% Scenario 3
     
    case 'high_TFP' %
     title_g= { '% Welfare changes High TFP & Cur. Policy in Scenario 3'};
    load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36.dat'] )
  newpolt  =trantf_muprop_sec2_tfp1_gz07_gn_grate36;
      
   case 'high_TFP_RR50'; %
       title_g= { '% Welfare changes RR50 & Cur. Policy in Scenario of High TFP'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36;   
      
   case 'high_TFP_MED10'
     title_g= { '% Welfare changes Med copay 10% up & Cur. Policy in Scenario of High TFP'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36_med10.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36_med10;         
      
    case 'high_TFP_Ltc10'
     title_g= { '% Welfare changes LTC copay 10% up & Cur. Policy in Scenario of High TFP'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36_ltc10.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36_ltc10; 
      
   case 'high_TFP_Med_Ltc30'
     title_g= { '% Welfare changes Med and LTC copay 10% up vs Current Policy'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_high_tfp.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_high_tfp;
     
     
 case 'high_TFP_RR50_Ltc10'
 title_g= { '%  Welfare change RR50% & LTC copay 10% down & Cur. Policy in Scenario 3 '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36_LTC.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36_LTC;
 
 case 'high_TFP_RR50_MED10'
 title_g= { '%  Welfare change RR50% & MED copay 10% down & Cur. Policy in Scenario 3 '};
 load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36_MED.dat'] )
 newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36_MED;
          
  case 'CurPolicy_high_TFP'; %
      title_g= { '% Welfare High TFP vs Current Policy'};
    load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36;
 
%   case 'TFP06'
% title_g= { '% Welfare changes TFP Growth Rate 0.6% vs 0.3% '};
%  load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate6.dat'] )
%  newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate6;
   
 %% Scenario 3
%    case 'high_TFP_RR50'
%  title_g= { '%  Welfare change RR50% & Cur. Policy in Scenario 3 '};
%  load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36.dat'] )
%  newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36;

 end

 
%%  current policy 
   load( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3.dat'] )
          curpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3;
% 
%  if Scenario ==1
%            load( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3.dat'] )
%           curpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3;
%   elseif   Scenario ==2
%          load( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire.dat'] )
%          curpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire;
%   elseif   Scenario ==3
%         load( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36.dat'] )
%           curpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36;
%   end
% 
%  case 'CurPolicy_no_pop'; %
%     title_g= { '% Welfare changes LTC copay 10% vs Current Policy'};
%     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_no_pop.dat'] )
%     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_no_pop;
%   
%     case 'immig'; %
%      title_g= { '% Welfare changes LTC copay 10% up vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_immig.dat'] )
%       newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_immig;
%     
%     case 'immig_RR50'; %
%      title_g= { '% Welfare changes LTC copay 10% up vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_immig.dat'] )
%       newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_immig;
%       
%       case 'immig_MED10'; %
%      title_g= { '% Welfare changes LTC copay 10% up vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_med3_grate3_immig.dat'] )
%       newpolt =trantf_muprop_sec2_tfp1_gz07_gn_med3_grate3_immig;
%       
%         case 'immig_MED_LTC30'; %
%      title_g= { '% Welfare changes LTC copay 10% up vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_immig.dat'] )
%       newpolt =trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3_immig;
%       
%          case 'immig_LTC10'; %
%      title_g= { '% Welfare changes LTC copay 10% up vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_ltc3_grate3_immig.dat'] )
%       newpolt =trantf_muprop_sec2_tfp1_gz07_gn_ltc3_grate3_immig;
%      


%% ==============================================

co_h = [ 0 0 1;      
             1 0 0];
%         0 0 1;
%         1 0 0;
%       0 0.75 0.75;
%       0.75 0 0.75;
%       0.75 0.75 0;
%       0.25 0.25 0.25];

set(groot,'defaultAxesColorOrder',co_h)

line_h= {'-','-.',':',':'};
set(groot,'defaultAxesLineStyleOrder',line_h)  


%%
switch agent
    case  'all'
          [WG,Ujf1,Ujc1,~,logCjs,logCfs,Util] = gains(curpolt,newpolt);
    case 'worker' 
          [WG,Ujf1,Ujc1,~,logCjs,logCfs] = gains_worker(curpolt,newpolt);
   case 'retiree'   
         [WG,Ujf1,Ujc1,~,logCjs,logCfs] = gains_retiree(curpolt,newpolt);
end

% [WG,Ujf1,Ujc1] = gains(curpolt,curpolt);
aa = 0;
WG = real(WG(aa+1:end,:));
% WG=WG2;
%%
%     plot_util;
%     plot_Cons;
%%
 
windowSize=7;
 WG_MA = movavg(WG(:,:),type,windowSize);  % Moving Avara

%% グラフ
aa = 20; % From 100 age out of 120 Age to -60 age

%   range = [-100.5,3]; ee=9;

fig = figure('File','Welfare changes','Position',[100,100,900,450]);
    plot(WG_MA(1+aa:200,:), 'LineWidth',2)
    grid on
    line([0,200],[0,0],'Color','Black')
    line([80,80],[-300,300],'Color','Black','LineStyle','--')
    line([100-65,100-65],[-300,300],'Color','Black','LineStyle','--')
    line([80+10,80+10],[-300,300],'Color','Red','LineStyle','-')
    xlabel('Birth-Year Cohort','Fontsize',16)
    txt = ['Current', ' retirees'] ; text(1,ee,txt, 'Fontsize',16);
    txt1 = ['Current', ' workers'] ; text(40,ee,txt1, 'Fontsize',16);
    txt2 = ['          Future', '  cohorts '] ; text(95,ee,txt2, 'Fontsize',16);
    axis([0,180-aa*1, range] )
    set(gca,'XTickLabel',[100, 80, 60, 40, 20, 0, -20, -40, -60 ],'Fontsize',16)
    %legend('Low','Medium','High','Top 1%','Location','Southwest')
    legend('male-regular','female-regular','male-contingent','female-contingent','Location','southeast','Box','off','Fontsize',18);
    %legend('boxoff')
%    title( [''; char(title_g) ],  'Fontsize',14)
    set(gca,'Fontsize',18)
%    fig.PaperPositionMode = 'manual';
    orient(fig,'landscape')
    saveas(gcf, ['./Fig/welfare_', char(newpolicy) ],'epsc')
    saveas(gcf,['./Fig/welfare_', char(newpolicy) ],'jpg')
    print(fig,['./Fig/welfare_', char(newpolicy) ],'-dpdf')
    save( ['./output/walfare_' char(newpolicy) '_'  char(agent) '.mat'], 'WG_MA' );

% plot_comparison_agents

%% table: Calculation of Walfare for Retierment, Current & Future Generations  

StartAge = 20;   % From 100 Age 
EndAge  = 200;  % Upto -60 Age

WG_Table ={};   Mean_WG_S1={};   Mean_WG_S2={} ;  

% 厚生・原数値で作表
% WG_Table{1} = WG(1+StartAge:EndAge,:);  %  All Generations
% WG_Table{2} = WG(1+StartAge:1+StartAge+(100-65),:);  % Retirement from 100 Age to 65 Age
% WG_Table{3} = WG(1+StartAge+(100-64):1+StartAge+(100-20),:);  % Current From 64 Age to 20 Age
% WG_Table{4} = WG(1+StartAge+(100-19):EndAge,:);  % Future from 19 Age to -60 age

% 厚生・移動平均の値で作表
WG_Table{1} = WG_MA(2+StartAge:EndAge,:);  %  All Generations
WG_Table{2} = WG_MA(2+StartAge:1+StartAge+(100-65),:);  % Retirement from 100 Age to 65 Age
WG_Table{3} = WG_MA(1+StartAge+(100-64):1+StartAge+(100-20),:);  % Current From 64 Age to 20 Age
WG_Table{4} = WG_MA(1+StartAge+(100-19):EndAge,:);  % Future from 19 Age to -60 age

varNames =   {  'Agents  '    , 'Total',  'Retire', 'Current','Future' };
Agent_name =  [ "male-regular";  "female-regular"; "male-contingent"; "female-contingent" ];        

if Scenario ==1
       disp( [' Make  Table of Senario 1 ',char(newpolicy)]   ) 
     
        WG_S1 = WG_Table;
        for  i = 1:4
          Mean_WG_S1{i} = mean(WG_S1{i});
       end   
        
       Tab=table(Agent_name ,Mean_WG_S1{1}', ...
                           Mean_WG_S1{2}', Mean_WG_S1{3}', Mean_WG_S1{4}',...
                           'VariableNames',varNames);
       
%       writetable(Tab,'Table-mean-Senario1.csv')
       writetable(Tab, ['Table-mean-Senario1_',char(newpolicy),'.csv'])
 
else  %Scenario ==2
      disp( [' Make  Table of Senario 2 ',char(newpolicy)] ) 
    
       WG_S2 = WG_Table;
       for  i = 1:4
          Mean_WG_S2{i} = mean(WG_S2{i});
       end   
        
       Tab=table(Agent_name ,Mean_WG_S2{1}', ...
                           Mean_WG_S2{2}', Mean_WG_S2{3}', Mean_WG_S2{4}',...
                           'VariableNames',varNames);
       
%       writetable(Tab,'Table-mean-Senario2.csv')
       writetable(Tab, ['Table-mean-Senario2_',char(newpolicy),'.csv'])

       
end


end