%
% PLOTBG    plot balanced growth results from bg*.f90. For details,
%           see the appendix (../appendix/appendix.pdf).
%

%           Ellen McGrattan, 12-19-10
%           Revised, ERM, 2-11-15 

clc
clear all;
close all;



%%
 load ../bgt_JP.dat
 bgtcmp= bgt_JP;
 
 aa =20;
 bb = 105;

% load ../trantf_jp_inp_out/bgt_JP_byage_sec2_tfp1_gz07.dat
% bgtcmp = bgt_JP_byage_sec2_tfp1_gz07;


%% パラメターの出力要変更
iplot     = 1;
alb       = bgtcmp(1);
aub       = bgtcmp(2);
alpha     = bgtcmp(3);
beta      = bgtcmp(4);
delt1     = bgtcmp(5);
deli1     = bgtcmp(6);
delt2     = bgtcmp(7);
deli2     = bgtcmp(8);
gam       = bgtcmp(9);
gn        = bgtcmp(10);
gz        = bgtcmp(11);
lam       = bgtcmp(12);
tauc      = bgtcmp(13);
taud      = bgtcmp(14);
taul      = bgtcmp(15);
taupc     = bgtcmp(16);
taupu     = bgtcmp(17);
tfp       = bgtcmp(18);
thett1    = bgtcmp(19);
theti1    = bgtcmp(20);
thett2    = bgtcmp(21);
theti2    = bgtcmp(22);

B         = bgtcmp(23);
G         = bgtcmp(24);
alph1     = 1-alpha;
KT1       = bgtcmp(25);
KT2       = bgtcmp(26);
KI1       = bgtcmp(27);
KI2       = bgtcmp(28);
L1        = bgtcmp(29);
L2        = bgtcmp(30);
C         = bgtcmp(31);
XT1       = bgtcmp(32);
XT2       = bgtcmp(33);
XI1       = bgtcmp(34);
XI2       = bgtcmp(35);
Y         = bgtcmp(36);
ltax      = bgtcmp(37);
wpsi      = bgtcmp(38);
rpsi      = bgtcmp(39);
irate     = bgtcmp(40);
ctran     = bgtcmp(41);
nage      = bgtcmp(42);
nprod     = bgtcmp(43);
ntypes    = nage*nprod;
n         = bgtcmp(44);
thetl1    = 1-thett1-theti1;
thetl2    = 1-thett2-theti2;
rt1       = (irate-1)/(1-taupc)+delt1;
rt2       = (irate-1)/(1-taupu)+delt2;
ri1       = irate-1+deli1;
ri2       = irate-1+deli2;
w         = thetl1*alpha*Y/L1;
agrid     = linspace(alb,aub,n)';
eps       = bgtcmp(45:44+ntypes);
ps        = bgtcmp(45+ntypes:44+2*ntypes);
xi        = bgtcmp(45+2*ntypes:44+3*ntypes);
mu        = bgtcmp(45+3*ntypes:44+4*ntypes);
jpi       = zeros(nage,nprod);
cj        = zeros(nage,nprod);
lj        = zeros(nage,nprod);
aj        = zeros(nage,nprod);
apj       = zeros(nage,nprod);
ltaxj     = zeros(nage,nprod);
wpsij     = zeros(nage,nprod);
rpsij     = zeros(nage,nprod);
mtrj      = zeros(nage,nprod);
opt       = zeros(ntypes,n);
jpi(:)    = bgtcmp(45+4*ntypes:44+5*ntypes);
cj(:)     = bgtcmp(45+5*ntypes:44+6*ntypes);
lj(:)     = bgtcmp(45+6*ntypes:44+7*ntypes);
aj(:)     = bgtcmp(45+7*ntypes:44+8*ntypes);
apj(:)    = bgtcmp(45+8*ntypes:44+9*ntypes);
ltaxj(:)  = bgtcmp(45+9*ntypes:44+10*ntypes);
wpsij(:)  = bgtcmp(45+10*ntypes:44+11*ntypes);
rpsij(:)  = bgtcmp(45+11*ntypes:44+12*ntypes);
mtrj(:)   = bgtcmp(45+12*ntypes:44+13*ntypes);
opt(:)    = bgtcmp(45+13*ntypes:44+ntypes*(n+13));
ntx       = bgtcmp(45+ntypes*(n+13));  
kk        = 45+ntypes*(n+13);
wtax      = zeros(ntx,3);
wtax(:)   = bgtcmp(kk+1:kk+3*ntx);
kk        = kk+3*ntx+1;
ntr       = bgtcmp(kk);
rtra      = zeros(ntr,3);
rtra(:)   = bgtcmp(kk+1:kk+3*ntr);
kk        = kk+3*ntr+1;
nbeq      = bgtcmp(kk);
beq       = bgtcmp(kk+1:kk+nbeq);
jbeq      = bgtcmp(kk+nbeq+1:kk+2*nbeq);
jinh      = bgtcmp(kk+2*nbeq+1:kk+3*nbeq);
epsmat    = zeros(nage,nprod);
psmat     = zeros(nage,nprod);
ximat     = zeros(nage,nprod);
mumat     = zeros(nage,nprod);
epsmat(:) = eps;
psmat(:)  = ps;
ximat(:)  = xi;
mumat(:)  = mu;
psij      = wpsij+rpsij;
%pens     = bgtcmp(86151:86554); % add 2020/05/30　gridによって場所が変わってしまうので修正。
%health   = bgtcmp(86555:86958); % add 2020/05/30　gridによって場所が変わってしまうので修正。
bgtcmp_rev  = flipud(bgtcmp);
%以下追加：
zeta     = flipud(bgtcmp_rev(1221:1224));
reprate  = flipud(bgtcmp_rev(1217:1220));
gam      = flipud(bgtcmp_rev(1213:1216));
ltcare   = flipud(bgtcmp_rev(809:1212)); 
pens     = flipud(bgtcmp_rev(405:808)); 
health   = flipud(bgtcmp_rev(1:404));


for l=1:nbeq
  xi(jbeq(l))  = xi(jbeq(l))-beq(l);
  xi(jinh(l))  = xi(jinh(l))+mu(jbeq(l))*beq(l)/mu(jinh(l));
end

for j=1:nprod;
  tem        = find(epsmat(:,j)<1e-10);
  aret(j)    = tem(1);
  Tranj(j)   = sum(ximat(:,j).*mumat(:,j)+psij(:,j).*mumat(:,j));
  Tranwrk(j) = sum((ctran+ximat(1:aret(j)-1,j)+psij(1:aret(j)-1,j)).* ...
                    mumat(1:aret(j)-1,j));
  Tranret(j) = sum((ctran+ximat(aret(j):nage,j)+psij(aret(j):nage,j)).* ...
                    mumat(aret(j):nage,j));
end;

Y1        = tfp*KT1^thett1*KI1^theti1*L1^thetl1;
Y2        = tfp*KT2^thett2*KI2^theti2*L2^thetl2;
Ychk      = 2*(Y1)^alpha*(Y2)^alph1;
p1        = alpha*Y/Y1;
p2        = alph1*Y/Y2;
cprof     = p1*Y1-w*L1-delt1*KT1-XI1;
uprof     = p2*Y2-w*L2-delt2*KT2-XI2;
D         = p1*Y1-w*L1-XT1-XI1-taupc*cprof;
KT1eoy    = (1+gz)*(1+gn)*KT1;
KT2eoy    = (1+gz)*(1+gn)*KT2;
KI1eoy    = (1+gz)*(1+gn)*KI1;
KI2eoy    = (1+gz)*(1+gn)*KI2;
Beoy      = (1+gz)*(1+gn)*B;
vT1       = (1-taud)*KT1eoy;
vI1       = (1-taud)*(1-taupc)*KI1eoy;
vT2       = KT2eoy;
vI2       = (1-taupu)*KI2eoy;
v         = vT1+vT2+vI1+vI2;

Rev       = taud*D+taupc*cprof+taupu*uprof+ltax+tauc*C;
GDP       = Y-XI1-XI2;
Int       = (irate-1)*B;
Beq       = G+ctran+sum(Tranj)-Rev-Beoy+irate*B;
GovR      = Beq+Rev+max(-Int,0);
GovE      = G+ctran+sum(Tranj)+max(Int,0);
Def       = GovE-GovR;
NWhh      = v+Beoy;
NWgov     = KT1eoy+KT2eoy+KI1eoy+KI2eoy-v-Beoy;

Cj        = cj.*mumat;
Lj        = lj.*mumat;
Aj        = apj.*mumat;

%% 効用関数
%for j=1:nprod
%  U(j)    = sum((beta.^[0:nage-1]').*cumprod([1;psmat(1:nage-1,j)]).* ...
%              (log(cj(:,j))+gam*log(1-lj(:,j))));
%  Disc(j) = sum((beta.^[0:nage-1]').*cumprod([1;psmat(1:nage-1,j)]));
%end;

 for j=1:nprod
%     U(j) = sum( gam(j).*((1.-lj(:,j)).^(1.-zeta(j))-1.)./(1.-zeta(j)));
     if ( zeta(j) == 1 )
       U(j)   = sum((beta.^[0:nage-16]').*cumprod([1;psmat(1:nage-16,j)]).* ...
                (log(cj(1:end-15,j)) + gam(j) .* log(1.-lj(1:end-15,j))));
     else
       U(j)   = sum((beta.^[0:nage-16]').*cumprod([1;psmat(1:nage-16,j)]).* ...
                (log(cj(1:end-15,j)) + gam(j) .* ((1.-lj(1:end-15,j)).^(1.-zeta(j))-1.)./(1.-zeta(j))));
     end
  Disc(j) = sum((beta.^[0:nage-16]').*cumprod([1;psmat(1:nage-16,j)]));
 end;

age = 19+[1:nage];


disp('Results of OLG Model with Intangible Capital')
disp('--------------------------------------------')
disp('  Parameters')
disp('  ----------')
disp(sprintf('    Asset grid lower bound               %8.3f ',alb))
disp(sprintf('    Asset grid upper bound               %8.3f ',aub))
disp(sprintf('    Corporate share of output            %8.3f ',alpha))
disp(sprintf('    Discount factor                      %8.3f ',beta))
disp(sprintf('    Tangible depreciation, sector 1      %8.3f ',delt1))
disp(sprintf('    Intangible depreciation, sector 1    %8.3f ',deli1))
disp(sprintf('    Tangible depreciation, sector 2      %8.3f ',delt2))
disp(sprintf('    Intangible depreciation, sector 2    %8.3f ',deli2))
disp(sprintf('    Weight on leisure in utility         %8.3f ',gam))
disp(sprintf('    Growth of population (%%)             %8.3f ',gn*100))
disp(sprintf('    Growth of technology (%%)             %8.3f ',gz*100))
disp(sprintf('    Tax rate on consumption              %8.3f ',tauc))
disp(sprintf('    Tax rate on distributions            %8.3f ',taud))
disp(sprintf('    Tax rate on labor                    %8.3f ',taul))
disp(sprintf('    Tax rate on corporate profits        %8.3f ',taupc))
disp(sprintf('    Tax rate on noncorp. profits         %8.3f ',taupu))
disp(sprintf('    Tangible capital share, sector 1     %8.3f ',thett1))
disp(sprintf('    Intangible capital share, sector 1   %8.3f ',theti1))
disp(sprintf('    Tangible capital share, sector 2     %8.3f ',thett2))
disp(sprintf('    Intangible capital share, sector 2   %8.3f ',theti2))
disp(sprintf('    TFP parameter                        %8.3f ',tfp))
disp(        '    Populations, by productivity level  ')
for j=1:nprod;
disp(sprintf('      Group %g                            %8.3f ',[j,sum(mumat(:,j))]))
end;
disp(sprintf('    Fraction working                     %8.3f ',sum(mu.*(eps>0))))
for j=1:nprod;
disp(sprintf('      Group %g                            %8.3f ',[j,sum(mumat(:,j).*(epsmat(:,j)>0))]))
end;
disp(sprintf('    Fraction retired                     %8.3f ',1-sum(mu.*(eps>0))))
for j=1:nprod;
disp(sprintf('      Group %g                            %8.3f ',[j,sum(mumat(:,j).*(epsmat(:,j)==0))]))
end;
disp(sprintf('    Age of first job                     %8g ',20))
disp(        '    Ages when retired                    ')
for j=1:nprod;
disp(sprintf('      Group %g                            %8g ',[j,aret(j)+19]))
end;
disp(sprintf('    Maximum age                          %8g ',nage+19))
disp(sprintf('    Degree of annuitization              %8g ',lam*100))
disp(' ')
disp('  Labor Tax Schedule')
disp('  ------------------')
disp(        '    Over incomes (y): ')
for i=1:ntx
  if (i<10);
    disp(sprintf('      y%g                                 %8.3f',[i,wtax(i,1)]))
  else
    disp(sprintf('      y%g                                %8.3f',[i,wtax(i,1)]))
  end;
end;
disp(' ')
disp(        '    Tax schedule is (a+b*y): ')
for i=1:ntx
  if (i<10);
    disp(sprintf('      a%g                                 %8.3f',[i,wtax(i,2)]))
  else
    disp(sprintf('      a%g                                %8.3f',[i,wtax(i,2)]))
  end;
end;
disp(' ')
for i=1:ntx
  if (i<10);
    disp(sprintf('      b%g                                 %8.3f',[i,wtax(i,3)]))
  else
    disp(sprintf('      b%g                                %8.3f',[i,wtax(i,3)]))
  end;
end;
disp(' ')
disp('  Retiree Transfer Schedule')
disp('  -------------------------')
disp(        '    Over incomes (y): ')
for i=1:ntr
  disp(sprintf('      y%g                                 %8.3f',[i,rtra(i,1)]))
end;
disp(' ')
disp(        '    Transfer schedule is (a+b*y): ')
for i=1:ntr
  disp(sprintf('      a%g                                 %8.3f',[i,rtra(i,2)]))
end;
disp(' ')
for i=1:ntr
  disp(sprintf('      b%g                                 %8.3f',[i,rtra(i,3)]))
end;
disp(' ')
disp('  Steady State Values')
disp('  -------------------')
disp(sprintf('    Interest rate (%%)                    %8.3f ',100*(irate-1)))
disp(sprintf('    Wage rate                            %8.3f ',w))
disp(sprintf('    HH transfer                          %8.3f ',ctran))
disp(' ')
disp('  National Accounts (rel. to GDP)  ')
disp('  -------------------------------  ')
disp('  * Gross National Product/Income ')
disp(' ')  
disp(sprintf('    Consumption (C)                      %8.3f ',C/GDP))
disp(sprintf('    Tangible investment (XT)             %8.3f ',(XT1+XT2)/GDP))
disp(sprintf('      Corporate (XT1)                    %8.3f ',XT1/GDP))
disp(sprintf('      Noncorporate (XT2)                 %8.3f ',XT2/GDP))
disp(sprintf('    Government spending (G)              %8.3f ',G/GDP))
disp(     '                                            =====')
disp(sprintf('    GDP (Y-XI)                           %8.3f ',1.))
disp(' ')
disp(sprintf('    Labor income (wL)                    %8.3f ',w*(L1+L2)/GDP))
disp(sprintf('      Corporate (wL1)                    %8.3f ',w*L1/GDP))
disp(sprintf('      Noncorporate (wL2)                 %8.3f ',w*L2/GDP))
disp(sprintf('    Capital income (Y-wL-XI)             %8.3f ',(cprof+delt1*KT1+uprof+delt2*KT2)/GDP))
disp(sprintf('      Corporate (p1 Y1-wL1-XI1)          %8.3f ',(cprof+delt1*KT1)/GDP))
disp(sprintf('        Taxes                            %8.3f ',taupc*cprof/GDP))
disp(sprintf('        Net dividends                    %8.3f ',D/GDP))
disp(sprintf('        Undistributed                    %8.3f ',(XT1-delt1*KT1)/GDP))
disp(sprintf('        Depreciation                     %8.3f ',delt1*KT1/GDP))

disp(sprintf('      Noncorporate (p2 Y2-wL2-XI2)       %8.3f ',(uprof+delt2*KT2)/GDP))
disp(sprintf('        Taxes                            %8.3f ',taupu*uprof/GDP))
disp(sprintf('        Distributions before tax         %8.3f ',(1-taupu)*uprof/GDP))
disp(sprintf('        Depreciation                     %8.3f ',delt2*KT2/GDP))
disp(     '                                            =====')
disp(sprintf('    GDI = GDP                            %8.3f ',1.))
disp(sprintf('     Corporate value added               %8.3f ',(w*L1+cprof+delt1*KT1)/GDP))
disp(sprintf('     Noncorporate value added            %8.3f ',(w*L2+uprof+delt2*KT2)/GDP))
disp(' ')
disp(sprintf('    Addenda: intangible (XI)             %8.3f ',(XI1+XI2)/GDP))
disp(' ')
disp('  * Government Expenditures/Receipts ')
disp(' ')  
disp(sprintf('    Consumption of G&S                   %8.3f ',G/GDP)) 
disp(sprintf('    Transfers                            %8.3f ',(ctran+sum(Tranj))/GDP)) 
disp(sprintf('      To retirees                        %8.3f ',sum(Tranret)/GDP))
for j=1:nprod
disp(sprintf('        Group %g                          %8.3f',[j,Tranret(j)/GDP]))
end
disp(sprintf('      To workers                         %8.3f ',sum(Tranwrk)/GDP))
for j=1:nprod
disp(sprintf('        Group %g                          %8.3f',[j,Tranwrk(j)/GDP]))
end
disp(sprintf('    Interest paid                        %8.3f ',max(Int/GDP,0))) 
disp(     '                                            =====')
disp(sprintf('    Expenditures                         %8.3f ',GovE/GDP))
disp(' ')
disp(sprintf('    Tax revenues                         %8.3f ',Rev/GDP))
disp(sprintf('      Consumption                        %8.3f ',tauc*C/GDP))
disp(sprintf('      Labor                              %8.3f ',ltax/GDP))
for j=1:nprod
disp(sprintf('        Group %g                          %8.3f ',[j,sum(ltaxj(:,j).*mumat(:,j))/GDP]))
end;
disp(sprintf('      Dividends                          %8.3f ',taud*D/GDP))
disp(sprintf('      C-Profits                          %8.3f ',taupc*cprof/GDP))
disp(sprintf('      NC-Profits                         %8.3f ',taupu*uprof/GDP))
disp(sprintf('    Interest received                    %8.3f ',max(-Int/GDP,0))) 
disp(sprintf('    Accidental bequests                  %8.3f ',Beq/GDP))
disp(     '                                            =====')
disp(sprintf('    Receipts                             %8.3f ',GovR/GDP))
disp(     '                                            =====')
disp(sprintf('    Deficit=Expend.-Receipts             %8.3f ',Def/GDP))

disp(' ')
disp('  Fixed Assets (rel. to GDP)')
disp('  --------------------------')
disp(sprintf('    Tangible capital (KT'')               %8.3f ',(KT1eoy+KT2eoy)/GDP))
disp(sprintf('      Corporate (KT1'')                   %8.3f ',KT1eoy/GDP))
disp(sprintf('      Noncorporate (KT2'')                %8.3f ',KT2eoy/GDP))
disp(sprintf('    Addenda: intangible (KI'')            %8.3f ',(KI1eoy+KI2eoy)/GDP))
disp(sprintf('      Corporate (KI1'')                   %8.3f ',KI1eoy/GDP))
disp(sprintf('      Noncorporate (KI2'')                %8.3f ',KI2eoy/GDP))
disp(     '                                            =====')
disp(sprintf('    Total                                %8.3f ',(KT1eoy+KT2eoy+KI1eoy+KI2eoy)/GDP))
disp(' ')
disp('  Flow of Funds (rel. to GDP)')
disp('  ---------------------------')
disp('  * Household Balance Sheet')
disp(' ')
disp(sprintf('    Tangible capital (vT)                %8.3f ',(vT1+vT2)/GDP))
disp(sprintf('      Corporate (vT1)                    %8.3f ',vT1/GDP))
disp(sprintf('      Noncorporate (vT2)                 %8.3f ',vT2/GDP))
disp(sprintf('    Government Debt (B'')                 %8.3f ',Beoy/GDP))
disp(sprintf('    Addenda: intangible (vI)             %8.3f ',(vI1+vI2)/GDP))
disp(sprintf('      Corporate (vI1)                    %8.3f ',vI1/GDP))
disp(sprintf('      Noncorporate (vI2)                 %8.3f ',vI2/GDP))
disp(     '                                            =====')
disp(sprintf('    Net Worth                            %8.3f ',NWhh/GDP))
disp(' ')
disp('  * Government Balance Sheet')
disp(' ')
disp(sprintf('    Tangible capital (KT''-vT)            %8.3f ',(KT1eoy+KT2eoy-vT1-vT2)/GDP))
disp(sprintf('    Government Debt (-B'')                %8.3f ',-Beoy/GDP))
disp(sprintf('    Addenda: intangible (KI''-vI)         %8.3f ',(KI1eoy+KI2eoy-vI1-vI2)/GDP))
disp(     '                                            =====')
disp(sprintf('    Net Worth                            %8.3f ',NWgov/GDP))
disp(' ')
disp('  Welfare')
disp('  -------')
disp('  * Individual consumption')
disp(sprintf('      Average across all types           %8.3f ',mean(mean(cj))))
disp(sprintf('      Weighted sum across all types      %8.3f ',sum(sum(mumat.*cj))))
for j=1:nprod
disp(sprintf('        Sum mu_j*c_j, group %g            %8.3f ',[j,sum(mumat(:,j).*cj(:,j))]))
end;
disp(' ')
disp('  * Individual labor supply ')
disp(sprintf('      Average across all types           %8.3f ',mean(mean(lj))))
disp(sprintf('      Weighted sum across all types      %8.3f ',sum(sum(mumat.*lj))))
for j=1:nprod
disp(sprintf('        Sum mu_j*l_j, group %g            %8.3f ',[j,sum(mumat(:,j).*lj(:,j))]))
end;
disp(' ')
disp('  * Individual effective labor supply ')
disp(sprintf('      Average across all types           %8.3f ',mean(mean(lj.*epsmat))))
disp(sprintf('      Weighted sum across all types      %8.3f ',sum(sum(mumat.*lj.*epsmat))))
for j=1:nprod
disp(sprintf('        Sum mu_j*l_j*eps_j, group %g      %8.3f ',[j,sum(mumat(:,j).*lj(:,j).*epsmat(:,j))]))
end;
disp(' ')
disp('  * Individual labor incomes rel. to GDP')
disp(sprintf('      Weighted sum across all types      %8.3f ',w/GDP*sum(sum(mumat.*epsmat.*lj))))
for j=1:nprod
disp(sprintf('        Sum w*mu_j*l_j*eps_j, group %g    %8.3f ',[j,w/GDP*sum(mumat(:,j).*epsmat(:,j).*lj(:,j))]))
end;
disp(' ')
disp('  * Values for computing welfare')
disp(' ')

%効用出力
for j=1:nprod
disp(sprintf('    Ind. lifetime utility, group %g       %8.3f ',[j,U(j)]))
disp(sprintf('    Discount = sum beta^j*prod(ps(j)))   %8.3f ',Disc(j)))
disp(' ')
end;

disp(' ')
disp('    % Welfare gain=(exp((U1-U2)/Discount)-1)*100 ')
disp(' ')
disp(' ')
disp(sprintf('  Detrended GDP:                         %8.3f ',GDP))


%
% Construct statistics by taxable income bin
%
bins  = [[0:.3:2.1]',[[.3:.3:2.1],inf]'];
nbins = length(bins);
yj    = irate*aj+w*epsmat.*lj;
yjavg = sum(sum(mumat.*yj));
yjrat = yj/yjavg;

yrhst = zeros(nbins,1);
trhst = zeros(nbins,1);

for t=1:nage;
  for j=1:nprod;
    for i=1:nbins
      if (yjrat(t,j)>=bins(i,1) & yjrat(t,j)<bins(i,2)); 
        yrhst(i,1) = yrhst(i,1)+mumat(t,j);
      end
    end;
  end;
end;

%% make figures %%
if iplot==1;

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Cross-sectional Consumption over the life cycle','Position', [300 300 700 500]);
           plot(age,Cj,'LineWidth',2)
           xlabel('Age')
           grid on
           title('Cross-sectional Consumption over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent')
           set(gca, 'Fontsize',14)
           saveas(gcf,'./Fig/Cross_sec_comp','epsc')
           
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Cross-sectional asset holdings over the life cycle','Position', [300 300 700 500]);
           plot(age,Aj,'LineWidth',2)
           xlabel('Age')
           grid on
           title('Cross-sectional asset holdings over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent')
           set(gca, 'Fontsize',12)
           saveas(gcf,'./Fig/Cross_sec_asset','epsc')
           
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Cross-sectional labor supply over the life cycle','Position', [300 300 700 500]);
	       plot(age,Lj,'LineWidth',2)
           xlabel('Age')
           grid on
           title('Cross-sectional labor supply over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent')
           set(gca, 'Fontsize',12)
           saveas(gcf,'./Fig/Cross_sec_lab','epsc')
           
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Individual consumption over the life cycle','Position', [300 300 700 500]);
           plot(age,cj,'LineWidth',2)
           xlabel('Age')
           grid on
           title('Individual consumption over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent','Location','northwest')
           set(gca, 'Fontsize',12)     
           saveas(gcf,'./Fig/ind_consump','epsc')
           
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Individual asset holdings over the life cycle','Position', [300 300 750 550]);
           plot(age,apj,'LineWidth',2)
           xlabel('Age')
           grid on
%           title('Individual asset holdings over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent','box','off','Fontsize',15.5,'location','northwest')
           set(gca, 'Fontsize',16)
           saveas(gcf,'./Fig/calib_ind_asset','epsc')
           saveas(gcf,'./Fig/calib_ind_asset','jpg')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Individual labor supply over the life cycle','Position', [300 300 700 500]);
figure(6), plot(age,lj,'LineWidth',2)
           xlabel('Age')
           grid on
%          title('Individual labor supply over the life cycle')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent','box','off','Fontsize',18)
           set(gca, 'Fontsize',16)
           saveas(gcf,'./Fig/calib_ind_lab','epsc')
           saveas(gcf,'./Fig/calib_ind_lab','jpg')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Pop distribution by type','Position', [300 300 700 500]);
	plot(age, mu(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, mu(102:202),'LineWidth',2)
    plot(age, mu(203:303),'LineWidth',2)
    plot(age, mu(304:404),'LineWidth',2)
    title('poplation distribution by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent')
    set(gca, 'Fontsize',12)
    saveas(gcf,'./Fig/type_pop','epsc')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Efficiency by type','Position', [300 300 700 500]);
	plot(age, eps(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, eps(102:202),'LineWidth',2)
    plot(age, eps(203:303),'LineWidth',2)
    plot(age, eps(304:404),'LineWidth',2)
%    title('efficiency by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent', 'Box','off', 'Fontsize',18)
    set(gca, 'Fontsize',16)
    saveas(gcf,'./Fig/calib_type_effic','epsc')
    saveas(gcf,'./Fig/calib_type_effic','jpg')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Med&LTC by type','Position', [300 300 700 500]);
    l1=plot(age, health(1:101),'LineWidth',2);
    xlabel('Age')
    grid on
    hold on
%     l1=plot(age, health(102:202),'LineWidth',2)
%     plot(age, health(203:303),'LineWidth',2)
%     plot(age, health(304:404),'LineWidth',2)
    l2=plot(age, ltcare(1:101),'LineWidth',2);
%     plot(age, ltcare(102:202),'LineWidth',2)
%     plot(age, ltcare(203:303),'LineWidth',2)
%     plot(age, ltcare(304:404),'LineWidth',2)
    title('Individual health & long term care by types')
    axis tight;    xlim([aa, bb ]);
    legend([l1,l2],{'health','long term care'},'Location','northwest')
    set(gca, 'Fontsize',12)
    saveas(gcf,'./Fig/type_hlth','epsc')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Pension by type','Position', [300 300 700 500]);
	plot(age, pens(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, pens(102:202),'LineWidth',2)
    plot(age, pens(203:303),'LineWidth',2)
    plot(age, pens(304:404),'LineWidth',2)
%    title('Individual pens by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent','Location','northwest','box','off','Fontsize',18)
    set(gca, 'Fontsize',16)
    saveas(gcf,'./Fig/calib_ind_type_pen','epsc')
    saveas(gcf,'./Fig/calib_ind_type_pen','jpg')
    
 %%  aggregation
      Pens = mu .* pens;
      Health = mu .* health;
      Ltcare = mu .* ltcare;
%%      
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Health aggregate','Position', [300 300 700 500]);
	plot(age, Health(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, Health(102:202),'LineWidth',2)
    plot(age, Health(203:303),'LineWidth',2)
    plot(age, Health(304:404),'LineWidth',2)
    title('Aggregate Health by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent','Location','northwest')
    set(gca, 'Fontsize',12)
    saveas(gcf,'./Fig/type_hlth_agg','epsc')
    
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','LTC aggregate','Position', [300 300 700 500]);
    plot(age, Ltcare(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, Ltcare(102:202),'LineWidth',2)
    plot(age, Ltcare(203:303),'LineWidth',2)
    plot(age, Ltcare(304:404),'LineWidth',2)
    title('Aggregate Long term care by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent','Location','northwest')
    set(gca, 'Fontsize',12)  
    saveas(gcf,'./Fig/type_ltc_agg','epsc')
    
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Pension aggregate','Position', [300 300 700 500]);
    plot(age, Pens(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, Pens(102:202),'LineWidth',2)
    plot(age, Pens(203:303),'LineWidth',2)
    plot(age, Pens(304:404),'LineWidth',2)
    title('Aggregate Pens by types')
    axis tight;    xlim([aa, bb ]);
    legend('male-regular','female-regular','male-contingent','female-contingent','Location','northwest')
    set(gca, 'Fontsize',12)
    saveas(gcf,'./Fig/type_pen_agg','epsc')

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Survival rate','Position', [300 300 700 500]);
    plot(age, ps(1:101),'LineWidth',2)
    xlabel('Age')
    grid on
    hold on
    plot(age, ps(102:202),'--','LineWidth',2)
 %   title('survival rate')
    axis tight;    xlim([aa, bb ]);
    legend('male','female','Location','northwest','box','off', 'Fontsize',18)
    set(gca, 'Fontsize',16)    
    saveas(gcf,'./Fig/calib_suv_rate','epsc')
    saveas(gcf,'./Fig/calib_suv_rate','jpg')
    
end;

% shrs = [sum(mumat)',sum(mumat.*epsmat.*lj)'/sum(sum(mumat.*epsmat.*lj)), ...
%          Tranwrk(:)/sum(Tranwrk(:)),Tranret(:)/sum(Tranret(:))]'*100;

%人口のタイプ別シェアを出す
shrs = [  sum(mumat)',sum(mumat.*epsmat.*lj)'/sum(sum(mumat.*epsmat.*lj)),...
          Tranwrk(:)/sum(Tranwrk(:)),Tranret(1,:)'/sum(Tranret(1,:))]'*100;  


disp(' ')
disp(' ')
disp(' ')
disp(' % Share     % Share    % Share of Transfers')
disp('Population Compensation  Workers   Retirees')
disp('============================================')
disp(sprintf('%9.0f %9.0f  %9.0f %9.0f \n',shrs))

%% extra figres

health_mat = [health(1:101), health(102:202), health(203:303), health(304:404)];
ltcare_mat = [ltcare(1:101), ltcare(102:202), ltcare(203:303), ltcare(304:404)];

eps_mat = [eps(1:101), eps(102:202), eps(203:303), eps(304:404)];
income_mat = w * lj .* eps_mat;

%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Labor income','Position', [300 300 700 500]);
plot(age,income_mat,'LineWidth',2)
           xlabel('Age')
           grid on
%           title('Individual of labor income')
           axis tight;    xlim([aa, bb ]);
           legend('male-regular','female-regular','male-contingent','female-contingent','Location','northeast','box','off','Fontsize',18)
           set(gca, 'Fontsize',16) 
           saveas(gcf,'./Fig/calib_ind_labor_income','epsc')
           saveas(gcf,'./Fig/calib_ind_labor_income','jpg')

mu_all = mumat(:,1) + mumat(:,2) + mumat(:,3) + mumat(:,4);


%%             
%Position: 描画可能領域の位置とサイズ:[left bottom width height]
figure('Name','Pop distribution','Position', [300 300 700 500]);
     plot(age, mu_all,'LineWidth',2)
     xlabel('Age')
     grid on
     title('Total population distribution')
     axis tight;    xlim([aa, bb ]);  
 %    legend('male-regular','female-regular','male-contingent','female-contingent')
     set(gca, 'Fontsize',12)
     saveas(gcf,'./Fig/pop_dist','epsc')