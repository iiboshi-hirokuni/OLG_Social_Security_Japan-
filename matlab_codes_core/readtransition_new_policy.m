%
% PLOTTRAN    plot transition results 
% The user should search and replace the name of the .dat file 
% by the new one throughout the document
% load trantfWPHGhighertfp.dat
% clc
% close all

% load ../file_out/trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3.dat
% trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_med3_ltc3_grate3;

% 
% load ../file_out/trantf_muprop_sec2_tfp1_gz07_gn_med3_grate3.dat
% trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_med3_grate3;

%  load ../file_out/trantf_muprop_sec2_tfp1_gz07_gn_grate3.dat
%  trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_grate3;
 
% load ../file_out/trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3.dat
%  trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_rr50_grate3;

% load ../file_out/trantf_muprop_sec2_tfp1_gz07_gn_med3_grate6.dat
%  trantfWPHGhighertfp =trantf_muprop_sec2_tfp1_gz07_gn_med3_grate6;


%%=======================================
alb       = trantfWPHGhighertfp(1);
aub       = trantfWPHGhighertfp(2);
alpha     = trantfWPHGhighertfp(3);
beta      = trantfWPHGhighertfp(4);
delt1     = trantfWPHGhighertfp(5);
deli1     = trantfWPHGhighertfp(6);
delt2     = trantfWPHGhighertfp(7);
deli2     = trantfWPHGhighertfp(8);
gam       = trantfWPHGhighertfp(9);
gz        = trantfWPHGhighertfp(10);
lam       = trantfWPHGhighertfp(11);
thett1    = trantfWPHGhighertfp(12);
theti1    = trantfWPHGhighertfp(13);
thett2    = trantfWPHGhighertfp(14);
theti2    = trantfWPHGhighertfp(15);
B0        = trantfWPHGhighertfp(16);
alph1     = 1-alpha;
thetl1    = 1-thett1-theti1;
thetl2    = 1-thett2-theti2;
idebt     = trantfWPHGhighertfp(17);
ispend    = trantfWPHGhighertfp(18);
nj        = trantfWPHGhighertfp(19);
ne        = trantfWPHGhighertfp(20);
nh        = nj*ne;
T         = trantfWPHGhighertfp(21);
n         = trantfWPHGhighertfp(22);
ntx       = trantfWPHGhighertfp(23);
ntr       = trantfWPHGhighertfp(24);
nbeq      = trantfWPHGhighertfp(25);
tbeq      = trantfWPHGhighertfp(26);
eps       = trantfWPHGhighertfp(26+1:26+nh);
i         = 26+nh;
debtH     = trantfWPHGhighertfp(i+1:i+T);
ps        = zeros(T,nh);
pens      = zeros(T,nh);
health    = zeros(T,nh);
xi        = zeros(T,nh);
i         = i+T;
ps(:)     = trantfWPHGhighertfp(i+1:i+T*nh);
xi(:)     = trantfWPHGhighertfp(i+1+T*nh:i+2*T*nh);
pens(:)   = trantfWPHGhighertfp(i+1+2*T*nh:i+3*T*nh);
health(:) = trantfWPHGhighertfp(i+1+3*T*nh:i+4*T*nh);
i         = i+4*T*nh;
hpii      = trantfWPHGhighertfp(i+1:i+nh);
popi      = trantfWPHGhighertfp(i+1+nh:i+2*nh);
i         = i+2*nh;
Bn        = trantfWPHGhighertfp(i+1:i+T);
G         = trantfWPHGhighertfp(i+1+T:i+2*T);
tauc      = trantfWPHGhighertfp(i+1+2*T:i+3*T);  NP.tauc= tauc;
taud      = trantfWPHGhighertfp(i+1+3*T:i+4*T);
taul      = trantfWPHGhighertfp(i+1+4*T:i+5*T);
taupc     = trantfWPHGhighertfp(i+1+5*T:i+6*T);
taupu     = trantfWPHGhighertfp(i+1+6*T:i+7*T);
tfp       = trantfWPHGhighertfp(i+1+7*T:i+8*T);
gnt       = trantfWPHGhighertfp(i+1+8*T:i+9*T);
i         = i+9*T;
wtaxt     = zeros(T,3*ntx);
rtrat     = zeros(T,3*ntr);
wtaxt(:)  = trantfWPHGhighertfp(i+1:i+3*T*ntx);
rtrat(:)  = trantfWPHGhighertfp(i+1+3*T*ntx:i+3*T*(ntx+ntr));
i         = i+3*T*(ntx+ntr);
KT1       = trantfWPHGhighertfp(i+1:i+T);        NP.KT1= KT1;
KT2       = trantfWPHGhighertfp(i+1+T:i+2*T);    NP.KT2= KT2;
KI1       = trantfWPHGhighertfp(i+1+2*T:i+3*T);  NP.KI1= KI1;
KI2       = trantfWPHGhighertfp(i+1+3*T:i+4*T);  NP.KI2= KI2;
L1        = trantfWPHGhighertfp(i+1+4*T:i+5*T);   NP.L1= L1;
L2        = trantfWPHGhighertfp(i+1+5*T:i+6*T);    NP.L2= L2;
C         = trantfWPHGhighertfp(i+1+6*T:i+7*T); NP.C= C;
Y         = trantfWPHGhighertfp(i+1+7*T:i+8*T);  NP.Y= Y;
XT1       = trantfWPHGhighertfp(i+1+8*T:i+9*T);
XT2       = trantfWPHGhighertfp(i+1+9*T:i+10*T);
XI1       = trantfWPHGhighertfp(i+1+10*T:i+11*T);
XI2       = trantfWPHGhighertfp(i+1+11*T:i+12*T);
cprof     = trantfWPHGhighertfp(i+1+12*T:i+13*T);
uprof     = trantfWPHGhighertfp(i+1+13*T:i+14*T);
D         = trantfWPHGhighertfp(i+1+14*T:i+15*T);
GDP       = trantfWPHGhighertfp(i+1+15*T:i+16*T);  NP.GDP= GDP;
Bnew      = trantfWPHGhighertfp(i+1+16*T:i+17*T);
Ltax      = trantfWPHGhighertfp(i+1+17*T:i+18*T);  NP.Ltax= Ltax;
Psi       = trantfWPHGhighertfp(i+1+18*T:i+19*T);
irate     = trantfWPHGhighertfp(i+1+19*T:i+20*T);  NP.irate= irate;
w         = trantfWPHGhighertfp(i+1+20*T:i+21*T);  NP.w= w;
Tran      = trantfWPHGhighertfp(i+1+21*T:i+22*T);
GDP2      = trantfWPHGhighertfp(i+1+22*T:i+23*T);
a         = zeros(nh,T);
c         = zeros(nh,T);
l         = zeros(nh,T);
at        = zeros(nh,T);
i         = i+23*T;

disp(['To plotgain =' num2str(i+1+T*nh)   ])
% i+1+T*nh  %--> plotgain

a(:)      = trantfWPHGhighertfp(i+1:i+T*nh);
c(:)      = trantfWPHGhighertfp(i+1+T*nh:i+2*T*nh);
l(:)      = trantfWPHGhighertfp(i+1+2*T*nh:i+3*T*nh);
at(:)     = trantfWPHGhighertfp(i+1+3*T*nh:i+4*T*nh);
beq       = zeros(tbeq,nbeq);
jbeq      = zeros(tbeq,nbeq);
jinh      = zeros(tbeq,nbeq);
i         = i+4*T*nh;
if (nbeq>0);
  for t=1:tbeq;
    for j=1:nbeq;
       beq(t,j)  = trantfWPHGhighertfp(i+1);
       jbeq(t,j) = trantfWPHGhighertfp(i+2);
       jinh(t,j) = trantfWPHGhighertfp(i+3);
       i         = i+3;
    end;
  end;
end;
i=i+1;
for j=1:ne;
    for t=1:T;
ltaxj(t,j)=trantfWPHGhighertfp(i);
i=i+1;
end;
end;

for j=1:ne;
for t=1:T;
wpsij(t,j)=trantfWPHGhighertfp(i);
i=i+1;
end;
end;

for j=1:ne;
for t=1:T;
rpsij(t,j)=trantfWPHGhighertfp(i);
i=i+1;
end;
end;
i
for t=1:T;
agedep(t)=100*(trantfWPHGhighertfp(i));
i=i+1;
end;


for j=1:nh;
for t=1:T;
income(t,j)=trantfWPHGhighertfp(i);
i=i+1;
end;
end;

% i
% 
% for j=1:nh;
% for t=1:T;
% incomeafter(t,j)=trantfWPHGhighertfpR(i);
% i=i+1;
% 
% end;
% end;



psij      = wpsij+rpsij;

epsmat    = zeros(nj,ne);
epsmat(:) = eps;
if idebt==1;
  Bn      = Bn.*[GDP(2:T);GDP(T)];
end;
B         = [B0;Bn(1:T-1)]; NP.B= B;
if ispend==1;
  G       = G.*GDP;  NP.G= G;
end;
L         = L1+L2;
taudl     = [taud(1);taud(1:T-1)];
taupcl    = [taupc(1);taupc(1:T-1)];
taupul    = [taupu(1);taupu(1:T-1)];
rT1       = ((1-taudl).*irate./(1-taud)-1)./(1-taupc)+delt1;
rT2       = (irate-1)./(1-taupu)+delt2;
rI1       = (1-taudl).*(1-taupcl).*irate./((1-taud).*(1-taupc))-1+deli1;
rI2       = (1-taupul).*irate./(1-taupu)-1+deli2;
agrid     = linspace(alb,aub,n)';

popt      = zeros(T,nh);
popt(1,:) = popi(:)';
for i=1:ne;
  for t=2:T
    h          = (i-1)*nj+1;
    popt(t,h)  = (1.+gnt(t-1))*popt(t-1,h);
  %   popt(t,h)  = popt(t-1,h);
 end;
  for t=2:T;
    for j=2:nj;
      h        = (i-1)*nj+j;
      popt(t,h)= popt(t-1,h-1)*ps(t-1,h-1);
    % popt(t,h)= popt(t-1,h);
    end
  end
end;



bxi       = xi;
if (nbeq>0);
%   for t=1:T;
%     for i=1:nbeq;
%       if (tbeq==1); 
%         jbeq(1,i)
%         bxi(t,jbeq(1,i)) = bxi(t,jbeq(1,i))-beq(1,i);
%         bxi(t,jinh(1,i)) = bxi(t,jinh(1,i))+popt(t,jbeq(1,i))*beq(1,i)/ ...
%                             popt(t,jinh(1,i));
%       else
%         bxi(t,jbeq(t,i)) = bxi(t,jbeq(t,i))-beq(t,i);
%         bxi(t,jinh(t,i)) = bxi(t,jinh(t,i))+popt(t,jbeq(t,i))*beq(t,i)/ ...
%                             popt(t,jinh(t,i));
%       end
%     end
%   end
end;

Tranj            = zeros(T,ne);
Healthexp        = zeros(T,ne);
Trannret         = zeros(T,ne);
Tranret          = zeros(T,ne);
for i=1:ne;
  h              = (i-1)*nj+1:i*nj;
  tem            = find(eps(h)<1e-10);
  Jr(i)          = tem(1);
  Tranj(:,i)     = sum((bxi(:,h).*popt(:,h))')'./sum(popt')'+...
                   sum((pens(:,h).*popt(:,h))')'./sum(popt')'+...
                   sum((health(:,h).*popt(:,h))')'./sum(popt')'+psij(:,i);
  Healthexp(:,i)     = sum((health(:,h).*popt(:,h))')'./sum(popt')';
  Trannret(:,i)  = Tran.*sum(popt(:,(i-1)*nj+1:(i-1)*nj+Jr(i)-1)')'./ ...
                   sum(popt')'+sum((bxi(:,(i-1)*nj+1:(i-1)*nj+Jr(i)-1).* ...
                   popt(:,(i-1)*nj+1:(i-1)*nj+Jr(i)-1))')'./sum(popt')'+ ...
                   wpsij(:,i);
  Tranret(:,i)   = Tran.*sum(popt(:,(i-1)*nj+Jr(i):i*nj)')'./sum(popt')'+ ...
                   sum((bxi(:,(i-1)*nj+Jr(i):i*nj).* ...
                   popt(:,(i-1)*nj+Jr(i):i*nj))')'./sum(popt')'+...
                   sum((pens(:,(i-1)*nj+Jr(i):i*nj).* ...
                   popt(:,(i-1)*nj+Jr(i):i*nj))')'./sum(popt')'+...
                   sum((health(:,(i-1)*nj+Jr(i):i*nj).* ...
                   popt(:,(i-1)*nj+Jr(i):i*nj))')'./sum(popt')'+rpsij(:,i);
  Nwrk(:,i)      = sum(popt(:,(i-1)*nj+1:(i-1)*nj+Jr(i)-1)')'./sum(popt')';
  Nret(:,i)      = sum(popt(:,(i-1)*nj+Jr(i):i*nj)')'./sum(popt')';
end;

KT1eoy    = XT1+(1-delt1)*KT1;
KT2eoy    = XT2+(1-delt2)*KT2;
KI1eoy    = XI1+(1-deli1)*KI1;
KI2eoy    = XI2+(1-deli2)*KI2;
Beoy      = Bnew+irate.*B;

vT1       = (1-taud).*KT1eoy;
vI1       = (1-taud).*(1-taupc).*KI1eoy;
vT2       = KT2eoy;
vI2       = (1-taupu).*KI2eoy;
v         = vT1+vT2+vI1+vI2;
Rev       = taud.*D+taupc.*cprof+taupu.*uprof+Ltax+tauc.*(C-sum(Healthexp')');
Int       = (1.01-1).*B;
Beq       = G+Tran+sum(Tranj')'-Rev-Bnew;
GovR      = Beq+Rev+max(-Int,0);
GovE      = G+Tran+sum(Tranj')'+max(Int,0);
Def       = GovE-GovR;
NWhh      = v+Beoy;
NWgov     = KT1eoy+KT2eoy+KI1eoy+KI1eoy-v-Beoy;
GDP3      = Y-XI1-XI2;

a         = a';
c         = c';

l         = l';
at        = at';
aggC      = sum((popt.*c)')'./sum(popt')'; NP.aggC= aggC;
aggL      = sum((popt.*l.*(ones(T,1)*eps(:)'))')'./sum(popt')'; NP.aggL= aggL;
popul     = sum(popt');
popul     = popul';
aggA      = sum((popt.*a)')'./sum(popt')'; NP.aggA= aggA;
aggAt     = sum((popt.*at)')'./sum(popt')';
tkap      = sum((popt.*at.*(1-ps))')'./sum(popt')';

%%
c(1,:) = 0.25;
l(1,:) = 0;
%%


Uj        = zeros(T,ne);
Disc      = zeros(T,ne);
Frac      = zeros(T,ne);
for i=1:ne;
  s       = int2str(i);
  eval(['Cj',s,'= zeros(T,nj);']);
  eval(['Aj',s,'= zeros(T,nj);']);
  eval(['Lj',s,'= zeros(T,nj);']);
  eval(['Atj',s,'= zeros(T,nj);']);
  for t=1:T-nj+1;
    eval(['Aj',s,'(t,:)  = diag(a(t:t+nj-1,(i-1)*nj+1:i*nj));'])
    eval(['Cj',s,'(t,:)  = diag(c(t:t+nj-1,(i-1)*nj+1:i*nj));'])
    eval(['Lj',s,'(t,:)  = diag(l(t:t+nj-1,(i-1)*nj+1:i*nj));'])
    eval(['Atj',s,'(t,:) = diag(at(t:t+nj-1,(i-1)*nj+1:i*nj));'])
    eval(['cj= log(Cj',s,'(t,:));']);
    eval(['lej= log(1-Lj',s,'(t,:));']);
    Uj(t,i)  = sum((beta.^[0:nj-1]').* ...
               cumprod([1;diag(ps(t:t+nj-2,(i-1)*nj+1:i*nj-1))]).* ...
              (cj'+gam*lej')); % (cj'+gam*lej')); 
    Disc(t,i)= sum((beta.^[0:nj-1]').* ...
               cumprod([1;diag(ps(t:t+nj-2,(i-1)*nj+1:i*nj-1))]));
  end
end;
for i=1:ne;
  for t=1:T;
    Frac(t,i)= sum(popt(t,(i-1)*nj+1:i*nj)'.*(epsmat(:,i)>0))/sum(popt(t,:));
  end
end;

for i=1:ne;
  s     = int2str(i);
  eval(['Aj',s,'(T-nj+2:T,:) = ones(nj-1,1)*Aj',s,'(T-nj+1,:);'])
  eval(['Cj',s,'(T-nj+2:T,:) = ones(nj-1,1)*Cj',s,'(T-nj+1,:);'])
  eval(['Lj',s,'(T-nj+2:T,:) = ones(nj-1,1)*Lj',s,'(T-nj+1,:);'])
  eval(['Atj',s,'(T-nj+2:T,:)= ones(nj-1,1)*Atj',s,'(T-nj+1,:);'])
  Uj(T-nj+2:T,i)   = ones(nj-1,1)*Uj(T-nj+1,i);
  Disc(T-nj+2:T,i) = ones(nj-1,1)*Disc(T-nj+1,i);
end;

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
disp(sprintf('    Growth of technology (%%)             %8.3f ',gz*100))
disp(sprintf('    Tangible capital share, sector 1     %8.3f ',thett1))
disp(sprintf('    Intangible capital share, sector 1   %8.3f ',theti1))
disp(sprintf('    Tangible capital share, sector 2     %8.3f ',thett2))
disp(sprintf('    Intangible capital share, sector 2   %8.3f ',theti2))
disp(sprintf('    Age of first job                     %8g ',21))
disp(sprintf('    Age when retired                         '))
for i=1:ne
disp(sprintf(['      Group ',int2str(i), ...
              '                            %8g '],Jr(i)+20))
end
disp(sprintf('    Maximum age                          %8g ',nj+20))
disp(sprintf('    Degree of annuitization              %8g ',lam*100))
disp(' ')

for i=1:ne;
  ij           = (i-1)*nj+1:i*nj;
  aggAj(:,i)   = sum((popt(:,ij).*a(:,ij))')'./sum(popt')';
  aggCj(:,i)   = sum((popt(:,ij).*c(:,ij))')'./sum(popt')';
  aggLj(:,i)   = sum((popt(:,ij).*l(:,ij).*(ones(T,1)*eps(ij)'))')'./sum(popt')';
end;
AVG = [B,G,tauc,taud,taul,taupc,taupu,tfp, ...
       100*(irate-1),w,Tran, ...
       C./GDP,(XT1+XT2)./GDP,XT1./GDP,XT2./GDP,G./GDP, ...
       (delt1*KT1+delt2*KT2)./GDP,delt1*KT1./GDP,delt2*KT2./GDP, ...
          w.*(L1+L2)./GDP,cprof./GDP,taupc.*cprof./GDP,D./GDP, ...
          (XT1-delt1*KT1)./GDP,uprof./GDP,taupu.*uprof./GDP, ...
          (1-taupu).*uprof./GDP,(XI1+XI2)./GDP,G./GDP,  ...
       (Tran+sum(Tranj')')./GDP,sum(Tranret')'./GDP,sum(Trannret')'./GDP, ...
          Tran./GDP,max(Int./GDP,0),GovE./GDP,Rev./GDP, ...
          tauc.*(C-sum(Healthexp')')./GDP,Ltax./GDP,taud.*D./GDP, ...
          taupc.*cprof./GDP,taupu.*uprof./GDP,max(-Int./GDP,0), ...
          Beq./GDP,GovR./GDP,Def./GDP, ...
       (KT1eoy+KT2eoy)./GDP,KT1eoy./GDP,KT2eoy./GDP,(KI1eoy+KI2eoy)./GDP, ...
          KI1eoy./GDP,KI2eoy./GDP,(KT1eoy+KT2eoy+KI1eoy+KI2eoy)./GDP, ...
       (vT1+vT2)./GDP,vT1./GDP,vT2./GDP,B./GDP,(vI1+vI2)./GDP,vI1./GDP, ...
          vI2./GDP,NWhh./GDP, ...
       (KT1eoy+KT2eoy-vT1-vT2)./GDP,-B./GDP, ... 
          (KI1eoy+KI2eoy-vI1-vI2)./GDP,NWgov./GDP, ...
       mean(Cj1')',C,mean(Lj1(:,1:Jr(1)-1)')',L,GDP,Uj,Disc,sum(Frac')', ...
          Frac,sum(Tranret')',Tranret,Tranret./(GDP*ones(1,ne)), ...
          Trannret./(GDP*ones(1,ne)),aggCj,aggLj];
 %AVG = [mean(AVG(1:10,:));
 AVG = [mean(AVG(2:10,:));
        mean(AVG(11:20,:));
        mean(AVG(21:30,:));
        mean(AVG(31:60,:));
        mean(AVG(61:89,:));
        mean(AVG(90:200,:));
        mean(AVG(201:240,:))];
AVG2 = ltaxj./(GDP*ones(1,ne));
%AVG2 =[mean(AVG2(1:10,:));
AVG2 =[mean(AVG2(2:10,:));
       mean(AVG2(11:20,:));
       mean(AVG2(21:30,:));
       mean(AVG2(31:60,:));
       mean(AVG2(61:89,:));
       mean(AVG2(90:200,:));
       mean(AVG2(201:240,:))];
   
disp('                                                   Period Averages  ')
disp('  ----------------                --------------------------------------------------          ') 
disp('  Exogenous Inputs                 2-10    11-20    21-30    31-60    61-89   90-200  BG Path')
disp('  ----------------                ------------------------------------------------------------')
disp(sprintf('    Government borrowing       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,1)))
disp(sprintf('    Government spending        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,2)))
        disp('    Tax rates')
disp(sprintf('      Consumption              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,3)))
disp(sprintf('      Distributions            %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,4)))
disp(sprintf('      Labor                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,5)))
disp(sprintf('      Corporate profits        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,6)))
disp(sprintf('      Noncorporate profts      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,7)))
disp(sprintf('    TFP parameters             %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,8)))
disp(sprintf('    Transfers to retirees      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,71+3*ne)))
for i=1:ne;
disp(sprintf(['      Group',int2str(i), ...
              '                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,71+3*ne+i)))
end
disp(sprintf('    Fraction working           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,70+2*ne)))
for i=1:ne;
disp(sprintf(['      Group',int2str(i), ...
              '                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,70+2*ne+i)))
end
disp(' ')
disp('  Endogenous Variables, averaged')
disp('  ------------------------------')
disp(sprintf('    Interest rate (%%)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,9)))
disp(sprintf('    Wage rate                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,10)))
disp(sprintf('    HH transfer                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,11)))
disp(' ')
disp('  National Accounts (rel. to GDP)  ')
disp('  -------------------------------  ')
disp('  * Gross National Product/Income ')
disp(' ')  
disp(sprintf('    Consumption (C)            %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,12)))
disp(sprintf('    Tangible investment (XT)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,13)))
disp(sprintf('      Corporate (XT1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,14)))
disp(sprintf('      Noncorporate (XT1)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,15)))
disp(sprintf('    Government Spending (G)    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,16)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    GDP (Y)                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',ones(1,7)))
disp(' ')
disp(sprintf('    Depreciation (delT KT)     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,17)))
disp(sprintf('      Corporate (delT KT1)     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,18)))
disp(sprintf('      Noncorporate (delT KT2)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,19)))
disp(sprintf('    Labor compensation (wL)    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,20)))
        disp('    Profits, sum pi Yi-wLi-depTi KTi-XIi')
disp(sprintf('      Corporate                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,21)))
disp(sprintf('        Taxes                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,22)))
disp(sprintf('        Net dividends          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,23)))
disp(sprintf('        Undistributed          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,24)))
disp(sprintf('      Noncorporate             %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,25)))
disp(sprintf('        Taxes                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,26)))
disp(sprintf('        Pre-tax distributions  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,27)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    GDI = GDP                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',ones(1,7)))
disp(' ')
disp(sprintf('    Addenda: intangible (XI)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,28)))

disp(' ')
disp('  * Government Expenditures/Receipts ')
disp(' ')  
disp(sprintf('    Consumption of G&S         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,29)))
disp(sprintf('    Transfers                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,30)))
disp(sprintf('      To retired               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,31)))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
             '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,71+4*ne+i)))
end;
disp(sprintf('      To non-retired           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,32)))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
             '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,71+5*ne+i)))
end;
disp(sprintf('    Interest paid              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,34)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    Expenditures               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,35)))
disp(' ')
disp(sprintf('    Tax revenues               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,36)))
disp(sprintf('      Consumption              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,37)))
disp(sprintf('      Labor                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,38)))
for i=1:ne;
disp(sprintf(['        Group ',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG2(:,i)))
end;
disp(sprintf('      Dividends                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,39)))
disp(sprintf('      C-Profits                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,40)))
disp(sprintf('      U-Profits                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,41)))
disp(sprintf('    Interest received          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,42)))
disp(sprintf('    Accidental bequests        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,43)))
disp(' ')
disp(sprintf('    Receipts                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,44)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    Deficit=Expend.-Receipts   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,45)))
disp(' ')
disp('  Fixed Assets (rel. to GDP)')
disp('  --------------------------')
disp(sprintf('    Tangible capital (KT'')     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,46)))
disp(sprintf('      Corporate (KT1'')         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,47)))
disp(sprintf('      Noncorporate (KT2'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,48)))
disp(sprintf('    Addenda: intangible (KI'')  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,49)))
disp(sprintf('      Corporate (KI1'')         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,50)))
disp(sprintf('      Noncorporate (KI2'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,51)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    Total                      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,52)))
disp(' ')
disp('  Flow of Funds (rel. to GDP)')
disp('  ---------------------------')
disp('  * Household Balance Sheet')
disp(' ')
disp(sprintf('    Tangible capital (vT)      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,53)))
disp(sprintf('      Corporate (vT1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,54)))
disp(sprintf('      Noncorporate (vT2)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,55)))
disp(sprintf('    Government Debt (B'')       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,56)))
disp(sprintf('    Addenda: intangible (vI)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,57)))
disp(sprintf('      Corporate (vI1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,58)))
disp(sprintf('      Noncorporate (vI2)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,59)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    Net Worth                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,60)))
disp(' ')
disp('  * Government Balance Sheet')
disp(' ')
disp(sprintf('    Tangible capital (KT''-vT)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,61)))
disp(sprintf('    Government Debt (-B'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,62)))
disp(sprintf('    Addenda: intang. (KI''-vI)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,63)))
disp(     '                                  =====    =====    =====    =====    =====    =====    =====')
disp(sprintf('    Net Worth                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,64)))
disp(' ')
disp('  Welfare')
disp('  -------')
disp('  * Individual consumption')
disp(sprintf('      Weighted sum mu_j*c_j    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,66)))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,71+6*ne+i)))
end;
disp('  * Individual labor supply ')
disp(sprintf('      Weighted sum mu_j*l_j*e_j%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,68)))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,71+7*ne+i)))
end;
disp('  * Values for computing welfare')
disp(sprintf('      Lifetime utility (U) '))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,69+i)))
end;
disp(sprintf('      Discount (sum beta^j ps(t,j)) '))
for i=1:ne
disp(sprintf(['        Group ',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f'],AVG(:,69+ne+i)))
end;
disp(' ')
disp('    % Welfare gain=(exp((U1-U2)/Discount)-1)*100 ')
disp(' ')
disp('  ---------------------------')
disp(' ')
disp('  NOTES: ')
disp(sprintf('    (1) level of GDP           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f',AVG(:,69)))
disp(' ')

save_text_file

%% mak graphs
% Plot_Graph
