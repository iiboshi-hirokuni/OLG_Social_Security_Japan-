% clear all
close all;

aa = 5;
bb = 206;
cc = 2020;         % start of year
dd = 2080;         % end  of year

% setting color 
co_h = [ 0 0 1;      
        1 0 0];
set(groot,'defaultAxesColorOrder',co_h)

line_h= {'-','--',':','--'};
set(groot,'defaultAxesLineStyleOrder',line_h)

% folder_name={'../file_out/grid1000/'} ;
folder_name={'../file_out/grid1000_New/'} ;

Scenario = 2;   % Baseline=1, retirement age ext=2, tfp high growth=3
%ここでは使わない% Policy = 4;  % Current=1, RR50=2, Med copay 10%up =3, LTC copay 10%up =4 

%%
if Scenario ==1;  % Baseline=1,
%title_g= { '% Welfare changes RR50% vs Current Policy'};
   leg1={'RR50%'};
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3.dat'] )
   newpolt=trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3; 

   trantfWPHGhighertfp = newpolt;
   readtransition_new_policy
%%
%title_g= { '% Welfare changes Med copay 10% up vs Current Policy'};
%   leg2={'RR50% & Med copay 10% down'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_MED.dat'] )
%   newpolt=trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_MED; 
   leg2={'Med copay 30%'}; %医療費自己負担30%へ
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_MED_up.dat'] )
   newpolt=trantf_muprop_sec2_tfp1_gz07_gn_grate3_MED_up ;
   trantfWPHGhighertfp = newpolt;
   readtransition_new_policy_2
%%
%title_g= { '% Welfare changes LTC OPR 30% vs Current Policy'};
%   leg3={'RR50% & LTC copay 10% down'};
%   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_LTC.dat'] )
%   newpolt=trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_LTC; 
   leg3={'LTC copay 30%'}; %介護費自己負担30%へ
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_LTC_up.dat'] )
   newpolt=trantf_muprop_sec2_tfp1_gz07_gn_grate3_LTC_up; 
   trantfWPHGhighertfp = newpolt;
   readtransition_new_policy_3

elseif   Scenario ==2;   %retirement age ext=2,
% title_g= { '% Welfare retire extension vs Current Policy'};   
     leg1={'RR50%'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_retire; 
     trantfWPHGhighertfp = newpolt;
     readtransition_new_policy
% 'Retire_ext_MED10'; %
     leg2={'Med copay 30%'};
     title_g= { '% Welfare retire extension vs Current Policy'};   
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_MED_up.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_MED_up;       
     trantfWPHGhighertfp = newpolt;
     readtransition_new_policy_2
% 'Retire_ext_LTC10'; %
      leg3={'LTC copay 30%'};
      title_g= { '% Welfare retire extension vs Current Policy'};   
      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_LTC_up.dat'] )
      newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire_LTC_up;  
      trantfWPHGhighertfp = newpolt;
      readtransition_new_policy_3

%%
% シナリオ3：High TFP は使わない
elseif   Scenario ==3;  % 
        leg1={'RR50%'};
% title_g= { '% Welfare changes LTC OPR 30% vs Current Policy'};
%      load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_high_tfp.dat'] )
%      newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3_high_tfp;       
        load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate36;       
     trantfWPHGhighertfp = newpolt;
     readtransition_new_policy
%   title_g= { '% Welfare changes Med and LTC OPR 30% vs Current Policy'};
     leg2={'Med copay 10% up'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36_med10.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36_med10;         
     trantfWPHGhighertfp = newpolt;
     readtransition_new_policy_2
% title_g= { '% Welfare changes Med and LTC OPR 30% vs Current Policy'};
     leg3={'LTC copay 10% up'};
     load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36_ltc10.dat'] )
     newpolt =trantf_muprop_sec2_tfp1_gz07_gn_grate36_ltc10;     
     trantfWPHGhighertfp = newpolt;
     readtransition_new_policy_3   
end

%% Current Policy

%    load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3.dat'] )
%    trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_grate3;
%    readtransition_cur_policy

if Scenario ==1;  % current=1,
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3.dat'] )
   trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_grate3;
   readtransition_cur_policy

elseif Scenario ==2;  % Retire ext =2,
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire.dat'] )
   trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_grate3_retire;
   readtransition_cur_policy

elseif Scenario ==3;  % high TFP = 3,
   load ( [ char(folder_name) 'trantf_muprop_sec2_tfp1_gz07_gn_grate36.dat'] )
   trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_grate36;
   readtransition_cur_policy
end


%% 

TT=[2015:2220];

% normalized  =  553/NP.GDP(6);
% normalized2 =  553/NP2.GDP(6);
% normalized3 =  553/NP3.GDP(6);
% normalized0 =  553/CUR.GDP(6);

%normalized  =  553/CUR.GDP(6);
%normalized2 =  553/CUR.GDP(6);
%normalized3 =  553/CUR.GDP(6);
%normalized0 =  553/CUR.GDP(6);

GDP2020 = 553;
% normalized  = GDP2020/CUR.GDP(6);
% normalized2 = GDP2020/CUR.GDP(6);
% normalized3 = GDP2020/CUR.GDP(6);
% normalized0 = GDP2020/CUR.GDP(6);

 normalized  =  GDP2020/NP.GDP(6);
 normalized2 =  GDP2020/NP2.GDP(6);
 normalized3 =  GDP2020/NP3.GDP(6);
 normalized0 =  GDP2020/CUR.GDP(6);


%% 実質GDP：Baselineの2020を1に基準化

fig = figure('Name','Y','File','Y')
    plot(TT(2:106),NP.GDP(2:106)*normalized /GDP2020, 'r-.','LineWidth',2);
    hold on
    grid on
    plot(TT(2:106),NP2.GDP(2:106)*normalized2 /GDP2020,'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.GDP(2:106)*normalized3 /GDP2020,'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.GDP(2:106)*normalized0 /GDP2020,'k-' ,'LineWidth',2);
    hold off
    ylim([0.55 1]); set(gca,'FontSize',18);
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Real GDP (2020 = 1)','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'}, 'location', 'best','box','off','FontSize',18);   
    set(gca,'FontSize',16)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_GDP_1','epsc')
%        saveas(gcf,'./Fig/tran_GDP_1', 'jpg')
        print(fig,'./Fig/tran_GDP_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_GDP_2','epsc')
%        saveas(gcf,'./Fig/tran_GDP_2', 'jpg')
        print(fig,'./Fig/tran_GDP_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_GDP_3','epsc')
%        saveas(gcf,'./Fig/tran_GDP_3', 'jpg')
        print(fig,'./Fig/tran_GDP_3','-dpdf')
    end

%% 人口推移（規模、億人）
poptAGG = zeros(240,1);
    for x = 1:240;
        poptAGG(x,1) = sum(popt(x, 1:404)) ;
    end
poptAGG_mf = zeros(240,1);
    for x = 1:240;
        poptAGG_mf(x,1) = sum(popt(x, 1:101)) ;
    end
poptAGG_ff = zeros(240,1);
    for x = 1:240;
        poptAGG_ff(x,1) = sum(popt(x, 102:202)) ;
    end
poptAGG_mp = zeros(240,1);
    for x = 1:240;
        poptAGG_mp(x,1) = sum(popt(x, 203:303)) ;
    end
poptAGG_fp = zeros(240,1);
    for x = 1:240;
        poptAGG_fp(x,1) = sum(popt(x, 304:404)) ;
    end

%億人単位の人口に変換
pop_2015 = 1.2314 ;

poptAGG = poptAGG .* pop_2015 ;
poptAGG_mf = poptAGG_mf .* pop_2015 ;
poptAGG_ff = poptAGG_ff .* pop_2015 ;
poptAGG_mp = poptAGG_mp .* pop_2015 ;
poptAGG_fp = poptAGG_fp .* pop_2015 ;

%% データセーブ
if Scenario == 1
   save_data_cur=[CUR.GDP, CUR.KT1, CUR.KT2, CUR.KI1, CUR.KI2, CUR.L1, CUR.L2, poptAGG];   
   writematrix(save_data_cur,'output_cur.xls')
  
  save_data_np=[NP.GDP, NP.KT1, NP.KT2, NP.KI1, NP.KI2, NP.L1, NP.L2, poptAGG];   
   writematrix(save_data_np,'output_np.xls')
   
    save_data_np2=[NP2.GDP, NP2.KT1, NP2.KT2, NP2.KI1, NP2.KI2, NP2.L1, NP2.L2, poptAGG];   
   writematrix(save_data_np2,'output_np2.xls')

     save_data_np3=[NP3.GDP, NP3.KT1, NP3.KT2, NP3.KI1, NP3.KI2, NP3.L1, NP3.L2, poptAGG];   
   writematrix(save_data_np3,'output_np3.xls')
else
     save_data_cur=[CUR.GDP, CUR.KT1, CUR.KT2, CUR.KI1, CUR.KI2, CUR.L1, CUR.L2, poptAGG];   
   writematrix(save_data_cur,'output2_cur.xls')
  
  save_data_np=[NP.GDP, NP.KT1, NP.KT2, NP.KI1, NP.KI2, NP.L1, NP.L2, poptAGG];   
   writematrix(save_data_np,'output2_np.xls')
   
    save_data_np2=[NP2.GDP, NP2.KT1, NP2.KT2, NP2.KI1, NP2.KI2, NP2.L1, NP2.L2, poptAGG];   
   writematrix(save_data_np2,'output2_np2.xls')

     save_data_np3=[NP3.GDP, NP3.KT1, NP3.KT2, NP3.KI1, NP3.KI2, NP3.L1, NP3.L2, poptAGG];   
   writematrix(save_data_np3,'output2_np3.xls')
end    

%% 一人当たり実質GDP：Baselineの2020を1に基準化
fig = figure('Name','per capita Y','File','per capita Y')
    plot(TT(2:106), (NP.GDP(2:106)*normalized ./poptAGG(2:106)) / (NP.GDP(6)*normalized / poptAGG(6)), 'r-.','LineWidth',2);
    hold on
    grid on
    plot(TT(2:106), (NP2.GDP(2:106)*normalized ./poptAGG(2:106)) / (NP.GDP(6)*normalized / poptAGG(6)),'b--' ,'LineWidth',2);
    plot(TT(2:106), (NP3.GDP(2:106)*normalized ./poptAGG(2:106)) / (NP.GDP(6)*normalized / poptAGG(6)),'g:' ,'LineWidth',2);
    plot(TT(2:106), (CUR.GDP(2:106)*normalized ./poptAGG(2:106)) / (NP.GDP(6)*normalized / poptAGG(6)),'k-' ,'LineWidth',2);
    hold off
    ylim([0.7 1.3]); set(gca,'FontSize',18);
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Per capita real GDP (2020 = 1)','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);   
    set(gca,'FontSize',16)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_GDPpercap_1','epsc')
        saveas(gcf,'./Fig/tran_GDPpercap_1', 'jpg')
        print(fig,'./Fig/tran_GDPpercap_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_GDPpercap_2','epsc')
       saveas(gcf,'./Fig/tran_GDPpercap_2', 'jpg')
        print(fig,'./Fig/tran_GDPpercap_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_GDPpercap_3','epsc')
       saveas(gcf,'./Fig/tran_GDPpercap_3', 'jpg')
        print(fig,'./Fig/tran_GDPpercap_3','-dpdf')
    end
    
%% 人口推移
fig = figure('Name','Popuration','File','Popuration')
    plot(TT(2:106), poptAGG(2:106), 'k','LineWidth',3);
    hold on
    grid on
    plot(TT(2:106),poptAGG_mf(2:106), 'b-','LineWidth',2);
    plot(TT(2:106),poptAGG_ff(2:106), 'r-','LineWidth',2);
    plot(TT(2:106),poptAGG_mp(2:106), 'b--','LineWidth',2);
    plot(TT(2:106),poptAGG_fp(2:106), 'r--','LineWidth',2);
        hold off
    ylim([0.0 1.3]); set(gca,'FontSize',18);
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Popuration','FontSize',18);
     legend('all', 'male-regular', 'female-regular', 'male-contingent', 'female-contingent', 'location', 'best','box','off','FontSize',18);   
    set(gca,'FontSize',16)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')   
    if Scenario == 1
        saveas(gcf,'./Fig/tran_poptAGG','epsc')
        saveas(gcf,'./Fig/tran_poptAGG', 'jpg')
        print(fig,'./Fig/tran_poptAGG','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_poptAGG','epsc')
        saveas(gcf,'./Fig/tran_poptAGG', 'jpg')
        print(fig,'./Fig/tran_poptAGG','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_poptAGG','epsc')
        saveas(gcf,'./Fig/tran_poptAGG', 'jpg')
        print(fig,'./Fig/tran_poptAGG','-dpdf')
    end
    
%% Aggregate GDPとPer capita GDPの成長率

fig = figure('Name','GDP growth rate','File','GDP growth rate')
plot(TT(3:106),log(NP.GDP(3:106)*normalized /GDP2020)-log(NP.GDP(2:105)*normalized /GDP2020), 'r-.','LineWidth', 2)
hold on
grid on
plot(TT(3:106),log(NP2.GDP(3:106)*normalized /GDP2020)-log(NP2.GDP(2:105)*normalized /GDP2020), 'b--.','LineWidth',2)
plot(TT(3:106),log(NP3.GDP(3:106)*normalized /GDP2020)-log(NP3.GDP(2:105)*normalized /GDP2020), 'g:','LineWidth',2)
plot(TT(3:106),log(CUR.GDP(3:106)*normalized /GDP2020)-log(CUR.GDP(2:105)*normalized /GDP2020), 'k-','LineWidth',2)
hold off
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('GDP growth rate','FontSize',18);
      legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);   
    set(gca,'FontSize',16)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_GDPgrowth_1','epsc')
        saveas(gcf,'./Fig/tran_GDPgrowth_1', 'jpg')
        print(fig,'./Fig/tran_GDPgrowth_1','-dpdf')
       elseif Scenario == 2
        saveas(gcf,'./Fig/tran_GDPgrowth_2','epsc')
        saveas(gcf,'./Fig/tran_GDPgrowth_2', 'jpg')
        print(fig,'./Fig/tran_GDPgrowth_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_GDPgrowth_3','epsc')
        saveas(gcf,'./Fig/tran_GDPgrowth_3', 'jpg')
        print(fig,'./Fig/tran_GDPgrowth_3','-dpdf')
    end
    
fig = figure('Name','Per capita GDP growth rate','File','Per capita GDP growth rate')
    plot(TT(3:106),log(NP.GDP(3:106)*normalized ./poptAGG(3:106)) / (NP.GDP(6)*normalized / poptAGG(6)) - log(NP.GDP(2:105)*normalized ./poptAGG(2:105)) / (NP.GDP(6)*normalized / poptAGG(6)), 'r-.','LineWidth', 2)
    hold on
    grid on
    plot(TT(3:106),log(NP2.GDP(3:106)*normalized ./poptAGG(3:106)) / (NP2.GDP(6)*normalized / poptAGG(6)) - log(NP2.GDP(2:105)*normalized ./poptAGG(2:105)) / (NP2.GDP(6)*normalized / poptAGG(6)), 'b--.','LineWidth',2)
    plot(TT(3:106),log(NP3.GDP(3:106)*normalized ./poptAGG(3:106)) / (NP3.GDP(6)*normalized / poptAGG(6)) - log(NP3.GDP(2:105)*normalized ./poptAGG(2:105)) / (NP3.GDP(6)*normalized / poptAGG(6)), 'g:','LineWidth',2)
    plot(TT(3:106),log(CUR.GDP(3:106)*normalized ./poptAGG(3:106)) / (CUR.GDP(6)*normalized / poptAGG(6)) - log(CUR.GDP(2:105)*normalized ./poptAGG(2:105)) / (CUR.GDP(6)*normalized / poptAGG(6)), 'k-','LineWidth',2)
    hold off
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Per capita GDP growth rate','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'}, 'location', 'best','box','off','FontSize',18);   
    set(gca,'FontSize',16)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_1','epsc')
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_1', 'jpg')
        print(fig,'./Fig/tran_GDPpercapita_growth_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_2','epsc')
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_2', 'jpg')
        print(fig,'./Fig/tran_GDPpercapita_growth_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_3','epsc')
        saveas(gcf,'./Fig/tran_GDPpercapita_growth_3', 'jpg')
        print(fig,'./Fig/tran_GDPpercapita_growth_3','-dpdf')
    end

%%
fig = figure('Name','Aggregate Asset to GDP','File','aggA')
    plot(TT(2:106), NP.aggA(2:106)*normalized ./NP.GDP(2:106)/normalized, 'r-.','LineWidth',2);
    hold on
    plot(TT(2:106),NP2.aggA(2:106)*normalized ./NP2.GDP(2:106)/normalized,'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.aggA(2:106)*normalized ./NP3.GDP(2:106)/normalized,'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.aggA(2:106)*normalized ./CUR.GDP(2:106)/normalized,'k-' ,'LineWidth',2);
    hold off
    grid on
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('aggregate Asset','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'}, 'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_aggA_1','epsc')
        saveas(gcf,'./Fig/tran_aggA_1', 'jpg')
        print(fig,'./Fig/tran_aggA_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_aggA_2','epsc')
        saveas(gcf,'./Fig/tran_aggA_2', 'jpg')
        print(fig,'./Fig/tran_aggA_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_aggA_3','epsc')
        saveas(gcf,'./Fig/tran_aggA_3', 'jpg')
        print(fig,'./Fig/tran_aggA_3','-dpdf')
    end
    
%%
fig = figure('Name','Gov Bond','File','Gov Bond')
    plot(TT(2:106),NP.B(2:106)*normalized,'r-.','LineWidth',2);
    hold on
    plot(TT(2:106),NP2.B(2:106)*normalized2,'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.B(2:106)*normalized3,'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.B(2:106)*normalized0,'k-' ,'LineWidth',2);
    hold off
    grid on
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Gov Bond','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'}, 'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_Gov_Bond_1','epsc')
        saveas(gcf,'./Fig/tran_Gov_Bond_1', 'jpg')
        print(fig,'./Fig/tran_Gov_Bond_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_Gov_Bond_2','epsc')
        saveas(gcf,'./Fig/tran_Gov_Bond_2', 'jpg')
        print(fig,'./Fig/tran_Gov_Bond_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_Gov_Bond_3','epsc')
        saveas(gcf,'./Fig/tran_Gov_Bond_3', 'jpg')
        print(fig,'./Fig/tran_Bond_3','-dpdf')
    end

    
%%
fig = figure('Name','aggregate Cons to GDP','File','aggC')
    plot(TT(2:106),NP.aggC(2:106)*normalized0./(NP.GDP(2:106)*normalized0), 'r-.','LineWidth',2);
    grid on
    hold on
    plot(TT(2:106),NP2.aggC(2:106)*normalized0./(NP2.GDP(2:106)*normalized0),'b--','LineWidth',2);
    plot(TT(2:106),NP3.aggC(2:106)*normalized0./(NP3.GDP(2:106)*normalized0),'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.aggC(2:106)*normalized0./(CUR.GDP(2:106)*normalized0),'k-' ,'LineWidth',2);
    hold off
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Consumption to GDP ratio','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_aggCons_1','epsc')
        saveas(gcf,'./Fig/tran_aggCons_1', 'jpg')
        print(fig,'./Fig/tran_aggCons_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_aggCons_2','epsc')
        saveas(gcf,'./Fig/tran_aggCons_2', 'jpg')
        print(fig,'./Fig/tran_aggCons_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_aggCons_3','epsc')
        saveas(gcf,'./Fig/tran_aggCons_3', 'jpg')
        print(fig,'./Fig/tran_aggCons_3','-dpdf')
    end
    
%%
% fig = figure('Name','aggregate Labor','File','aggL')
%     plot(TT(2:106),NP.aggL(2:106), 'LineWidth',2);
%      hold on
%     plot(TT(2:106),CUR.aggL(2:106), 'LineWidth',2);
%     hold off
%     grid on
%    xlim([cc dd]); set(gca,'FontSize',18);
%     title('aggregate Labor supply','FontSize',18);
%     legend('New Policy','Current Policy', 'location', 'best');
%     set(gca,'FontSize',18)
%     saveas(gcf,'./Fig/tran_aggLabor','epsc')
%%    fig.PaperPositionMode = 'manual'; 
%    orient(fig,'landscape')


%%
fig = figure('Name','Capital','File','K')
    plot(TT(2:106),[NP.KT1(2:106)+NP.KT2(2:106)+NP.KI1(2:106)+NP.KI2(2:106)]*normalized0./(NP.GDP(2:106)*normalized0), 'r-.','LineWidth',2);
    hold on
    plot(TT(2:106),[NP2.KT1(2:106)+NP2.KT2(2:106)+NP2.KI1(2:106)+NP2.KI2(2:106)]*normalized0 ./(NP2.GDP(2:106)*normalized0),'b--', 'LineWidth',2);
    plot(TT(2:106),[NP3.KT1(2:106)+NP3.KT2(2:106)+NP3.KI1(2:106)+NP3.KI2(2:106)]*normalized0 ./(NP3.GDP(2:106)*normalized0),'g:', 'LineWidth',2);
    plot(TT(2:106),[CUR.KT1(2:106)+CUR.KT2(2:106)+CUR.KI1(2:106)+CUR.KI2(2:106)]*normalized0 ./(CUR.GDP(2:106)*normalized0),'k-', 'LineWidth',2);
    hold off
    grid on
    hold on
%     plot(TT(2:106),NP.KT2(2:106), 'LineWidth',2);
%     plot(TT(2:106),KI1(2:106), 'LineWidth',2);
%     plot(TT(2:106),KI2(2:106), 'LineWidth',2);
    hold off
    xlim([cc dd]); set(gca,'FontSize',18);
%     legend({'KT1', 'KT2', 'KI1', 'KI2'},'Location','best','FontSize',10)
%    title('Capital to GDP ratio: Tangible & Intangible','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_Capital_1','epsc')
        saveas(gcf,'./Fig/tran_Capital_1', 'jpg')
        print(fig,'./Fig/tran_Capital_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_Capital_2','epsc')
        saveas(gcf,'./Fig/tran_Capital_2', 'jpg')
        print(fig,'./Fig/tran_Capital_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_Capital_3','epsc')
        saveas(gcf,'./Fig/tran_Capital_3', 'jpg')
        print(fig,'./Fig/tran_Capital_3','-dpdf')
    end

    
%%    
fig = figure('Name','Interest rate','File','irate')
    plot(TT(2:106),NP.irate(2:106), 'r-.','LineWidth',2);
    hold on
    plot(TT(2:106),NP2.irate(2:106),'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.irate(2:106),'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.irate(2:106),'k-' ,'LineWidth',2);
    hold off
    grid on
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Interest rate','FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_irate_1','epsc')
        saveas(gcf,'./Fig/tran_irate_1', 'jpg')
        print(fig,'./Fig/tran_irate_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_irate_2','epsc')
        saveas(gcf,'./Fig/tran_irate_2', 'jpg')
        print(fig,'./Fig/tran_irate_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_irate_3','epsc')
        saveas(gcf,'./Fig/tran_irate_3', 'jpg')
        print(fig,'./Fig/tran_irate_3','-dpdf')
    end



 %%   
fig = figure('Name','Wage rate','File','Wage rate')
    plot(TT(2:106),NP.w(2:106), 'r-.','LineWidth',2);
    hold on
    plot(TT(2:106),NP2.w(2:106),'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.w(2:106),'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.w(2:106),'k-' ,'LineWidth',2);
    
    hold off
    grid on
    xlim([cc dd]); set(gca,'FontSize',18);
%    title('Wage rate','FontSize',18); 
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off','FontSize',18);    
    set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_wage_1','epsc')
        saveas(gcf,'./Fig/tran_wage_1', 'jpg')
        print(fig,'./Fig/tran_wage_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_wage_2','epsc')
        saveas(gcf,'./Fig/tran_wage_2', 'jpg')
        print(fig,'./Fig/tran_wage_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_wage_3','epsc')
        saveas(gcf,'./Fig/tran_wage_3', 'jpg')
        print(fig,'./Fig/tran_wage_3','-dpdf')
    end

    
    
%% Government
%  figure('Name','SNA','File','SNA') 
%      plot(TT(2:106), C(2:106)./GDP(2:106),'b-','LineWidth',2);
%      grid on
%      hold on
%      plot(TT(2:106),(XT1(2:106) + XT2(2:106))./GDP(2:106),'r-','LineWidth',2);
%      plot(TT(2:106),XT1(2:106)./GDP(2:106),'--','LineWidth',2);
%      plot(TT(2:106),XT2(2:106)./GDP(2:106),'-.','LineWidth',2);
%      plot(TT(2:106),(XI1(2:106) + XI2(2:106))./GDP(2:106),'g-','LineWidth',2);
%      plot(TT(2:106),XI2(2:106)./GDP(2:106),'--.','LineWidth',2);
%      plot(TT(2:106),XI2(2:106)./GDP(2:106),'-.','LineWidth',2);
%      plot(TT(2:106),G(2:106)./GDP(2:106),'g:','LineWidth',2);
%      hold off
%     xlim([cc dd]); set(gca,'FontSize',18);
%      ylim([ -0.01  0.7])
%      title('National Account','FontSize',18)
%      legend({'C','XT','XT1','XT2','XI(参考)','XI1(参考)','XI2(参考)','G' },'Location','best', 'FontSize',10)
%      set(gca,'FontSize',18)
%     saveas(gcf,'./Fig/tran_SNA','epsc')

% %%
fig= figure('Name','Gov','File','Gov') 
%     plot(TT(2:106),NP.B(2:106), 'LineWidth',2);
%     grid on
%     hold on
%     plot(TT(2:106),NP.G(2:106), 'LineWidth',2);
     plot(TT(2:106),CUR.B(2:106), 'LineWidth',2);
     grid on
     hold on
     plot(TT(2:106),CUR.G(2:106), 'LineWidth',2);    
     hold off
     xlim([cc dd]); set(gca,'FontSize',18);
     legend('Bond','G.Spending','location', 'best','box','off','FontSize',18);   
%     title('Gov. Borrowing and Spending','FontSize',18)
     set(gca,'FontSize',18)
%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_GovSpend_1','epsc')
        saveas(gcf,'./Fig/tran_GovSpend_1', 'jpg')
        print(fig,'./Fig/tran_GovSpend_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_GovSpend_2','epsc')
        saveas(gcf,'./Fig/tran_GovSpend_2', 'jpg')
        print(fig,'./Fig/tran_GovSpend_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_GovSpend_3','epsc')
        saveas(gcf,'./Fig/tran_GovSpend_3', 'jpg')
        print(fig,'./Fig/tran_GovSpend_3','-dpdf')
    end
     
%%
% figure('Name','Tran','File','Tran') 
%     plot(TT(2:106),Tranj(2:106,1)+Tranj(2:106,2)+Tranj(2:106,3)+Tranj(2:106,4), ...
%     'LineWidth',2);
%     grid on
%     hold on
%     plot(TT(2:106),Trannret(2:106,1)+Trannret(2:106,2)+Trannret(2:106,3)+Trannret(2:106,4), ...
%     'LineWidth',2);
%     plot(TT(2:106),Tranret(2:106,1)+Tranret(2:106,2)+Tranret(2:106,3)+Tranret(2:106,4), ...
%     'LineWidth',2);
%    xlim([cc dd]); set(gca,'FontSize',18);
%     legend({'Tranj','Trannret','Tranret'},'Location','best','FontSize',10)     
%     title('Tran','FontSize',18)
% %     legend('New Policy','Current Policy');
%     set(gca,'FontSize',18)
%     %saveas(gcf,'./Fig/tran_aggA','epsc')
    
%%  Tax
fig = figure('Name','Tax rate','File','tax') 
    plot(TT(2:106),NP.tauc(2:106), 'r-.','LineWidth',2);
    grid on
    hold on
    plot(TT(2:106),NP2.tauc(2:106),'b--' ,'LineWidth',2);
    plot(TT(2:106),NP3.tauc(2:106),'g:' ,'LineWidth',2);
    plot(TT(2:106),CUR.tauc(2:106),'k-' ,'LineWidth',2);
    hold off
    ylim([0 0.5]); set(gca,'FontSize',18);
    xlim([cc dd]); set(gca,'FontSize',18);
    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off', 'FontSize',18);    
%    title('Consumption tax','FontSize',18)
    set(gca,'FontSize',16)
%%    fig.PaperPositionMode = 'manual'; 
    orient(fig,'landscape')
    if Scenario == 1
        saveas(gcf,'./Fig/tran_tax_1','epsc')
        saveas(gcf,'./Fig/tran_tax_1', 'jpg')
        print(fig,'./Fig/tran_tax_1','-dpdf')
        elseif Scenario == 2
        saveas(gcf,'./Fig/tran_tax_2','epsc')
        saveas(gcf,'./Fig/tran_tax_2', 'jpg')
        print(fig,'./Fig/tran_tax_2','-dpdf')
        elseif Scenario == 3
        saveas(gcf,'./Fig/tran_tax_3','epsc')
        saveas(gcf,'./Fig/tran_tax_3', 'jpg')
        print(fig,'./Fig/tran_tax_3','-dpdf')
    end


% 
% 
% % readtransition
% 
% readtransition_new_policy

%%  population growth

type = 'linear';
windowSize=8;
 MA_popgt = movavg(gnt(:,:),type,windowSize)*100;  % Moving Avara

fig = figure('Name','poplation groth rate','File','pop') 
%    plot(TT(2:106),gnt(2:106), 'r-.','LineWidth',2);
%    plot(TT(1:206),gnt(1:206), 'b-','LineWidth',2);
    plot(TT(2:206),MA_popgt(2:206), 'b-','LineWidth',2);
    grid on
    ylim([-3. 3]); set(gca,'FontSize',18);
%    xlim([cc dd]); set(gca,'FontSize',18);
    xlim([2015 2110]); set(gca,'FontSize',18);
    xlabel('year')
%    legend({ char(leg1),char(leg2),char(leg3), 'Current'},'location', 'best','box','off', 'FontSize',18);    
%    title('Consumption tax','FontSize',18)
    set(gca,'FontSize',16)
    orient(fig,'landscape')
    saveas(gcf,'./Fig/calib_pop_growth','epsc')
    saveas(gcf,'./Fig/calib_pop_growth', 'jpg')
    print(fig,'./Fig/calib_pop_growth','-dpdf')