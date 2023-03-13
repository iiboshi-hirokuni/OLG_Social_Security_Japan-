%% 生存確率かけた期待値の医療費と介護費のグラフを作成します。
% まず、「plotbg_type.m」でデータ作成してから実行する。

%要介護率のデータを読み込みます。
%LTCrate = readmatrix('LTCrate.csv');
LTCrate = readmatrix('LTCrate_linear.csv'); %線形補間済み
LTCrate = LTCrate(:,3:4);
%LTCrate(1,:) = [] % 1行目のcol nameを削除

%%

aa=21; bb=100;

% 期待医療費＝年齢別医療費×生存確率、期待介護費＝年齢別介護費×生存確率×要介護率
figure('Name','Med_LTC_by_gender','File','Med_LTC_by_gender','Position', [100 100 1000 600])
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
subplot(1,2,1)
   l1=plot(age, 10*health(1:101), 'k-',  'LineWidth',1.5);
   grid on
   hold on
     l2=plot(age, 10*health(1:101).*ps(1:101),'b-','LineWidth',2.5);
     l3=plot(age, 10*health(1:101).*ps(102:202),'r','LineWidth',2.5);
   hold off   
    xlabel('Age');  xlim([aa, bb ]);
    ylabel('Million Yen')
    ylim([-0 3.3]);
%   title('Data and expected eedical expenditure, per capita')
    legend({'data:recipients of Med','male:agent','female:agent'},'Location','best','Box','off', 'Fontsize',19)
    set(gca, 'Fontsize',16);
    
 subplot(1,2,2)   
   l2=plot(age,10*ltcare(1:101),'k-','LineWidth',1.5);
   grid on
   hold on
     l2=plot(age, 10*ltcare(1:101).*ps(1:101).*LTCrate(:,1),'b-','LineWidth',2.5);
     l3=plot(age, 10*ltcare(1:101).*ps(102:202).*LTCrate(:,2),'r','LineWidth',2.5);
   hold off   
 %    title('Data and expected LTC expenditure, per capita')
    axis tight;    xlim([aa, bb ]);
     xlabel('Age');  ylabel('Million Yen');
    ylim([-0 3.3]);
     legend({'data:recipients of LTC','male:agent','female:agent'},'Location','best','Box','off', 'Fontsize',19)
    set(gca, 'Fontsize',16);
    saveas(gcf,'./Fig/calib_type_exp-medltc_expend','epsc')
    saveas(gcf,'./Fig/calib_type_exp-medltc_expend','jpg')
 

%% 1人あたり医療・介護費用：岩本福井グラフ
%エクセルからデータを読み込む
data_medltcexp = xlsread('medl-ltc-IwaFuku.xlsx','sheet1','A1:C22');
X = categorical({'0-4' '5-9' '10-14' '15-19' '20-24' '25-29' '30-34' '35-39' '40-44' '45-49' '50-54' '55-59' '60-64' '65-69' '70-74' '75-79' '80-84' '85-89' '90-94' '95-99' '100-'});
X = reordercats(X,{'0-4' '5-9' '10-14' '15-19' '20-24' '25-29' '30-34' '35-39' '40-44' '45-49' '50-54' '55-59' '60-64' '65-69' '70-74' '75-79' '80-84' '85-89' '90-94' '95-99' '100-'});

figure('Name','Per capita medical and LTC expenditure')
bar(X, [data_medltcexp(:,1)/1000 data_medltcexp(:,2)/1000],1.0, 'grouped', 'LineStyle', 'none' )
%hold on
%bar(X, data_medltcexp(:,2)/1000, 'grouped', 'LineStyle', 'none', 'FaceColor', 'b' )
% title('Per capita medical and LIC expenditure','Fontsize',14)
 legend({'Med','LTC'},'Box','off', 'Location','best','Fontsize',14)
 set(gca,'Fontsize',12,'FontName','Arial');
 grid on
  ylabel('thousand yen');
  xlabel('age bracket');
 saveas(gcf,'./Fig/med-ltc_exp','epsc')
 saveas(gcf,'./Fig/med-ltc_exp','jpg')
