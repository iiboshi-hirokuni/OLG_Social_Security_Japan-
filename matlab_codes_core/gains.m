function [WG,Ufs,Ucs,Disc,logCjs,logCfs,Util] = gains(curpolt,newpolt);
% GAINS  compute gains by cohort for the current policy (curpolt)
%        and the new policy (newpolt) for the transition.
%
%        Ellen McGrattan, 9-23-13

windowSize = 7; % à⁄ìÆïΩãœÇÃä˙ä‘

init = 5;

%  labor supply
retire_fg =1;  % yes=1 or no=0  for exclusion of labour disutility  of  Age 65 to 70
gam           = 10; %curpolt(9);
%  gam           = 0;
zeta           =  [ 0.05; 0.03; 0.05; 0.03];

beta          = curpolt(4);
nj            = curpolt(19);
ne            = curpolt(20);
T             = curpolt(21);
nh            = ne*nj;
ps            = zeros(T,nh);
ps(:)         = curpolt(671:670+T*nh);
ntx           = curpolt(23);
ntr           = curpolt(24);
i1 = 504759+960+1;
% i1            = 503318;
cc            = zeros(nh,T);
lc            = zeros(nh,T);
% cc(:)         = curpolt(i1+1:i1+nh*T);
% lc(:)         = curpolt(i1+1+nh*T:i1+2*nh*T);

%% add  2021/07/25
i         = 26+nh;
i         = i+T;
i         = i+4*T*nh;
i         = i+2*nh;
i         = i+9*T;
i         = i+3*T*(ntx+ntr);
i         = i+23*T;
 cc(:)      = curpolt(i+1+T*nh:i+2*T*nh);
 lc(:)      = curpolt(i+1+2*T*nh:i+3*T*nh);
 
%%
 
%  cc =   cc(:,init+1:T) ; 
%  lc =   lc(:,init+1:T) ;
% size(cc)

cc            = cc';
lc             = lc';
% 
cc = cc(2:end,:);
lc = lc(2:end,:);

if  retire_fg==1
lc(:,46:50)=0.001*ones(239,5);  % for male full-time
lc(:,46+101:50+101)=0.001*ones(239,5);
lc(:,46+2*101:50+2*101)=0.001*ones(239,5);
lc(:,46+3*101:50+3*101)=0.001*ones(239,5);
end
%   cc(1,:) = 0.25;
%   lc(1,:) = 0;


save('data_cc_lc.mat','cc','lc')

%%
  T = T-1;
%%

for j=2:nj
  Discl(j,1)  = sum((beta.^[j-1:nj-1]').* ...
                cumprod([diag(ps(1:nj-j+1,j-1:nj-1))]));
end;
for t=1:T-nj+1;
  Disc(t,1)   = sum((beta.^[0:nj-1]').* ...
                cumprod([1;diag(ps(t:t+nj-1,1:nj-1))]));
end;
Disc(T-nj+2:T)= ones(nj-1,1)*Disc(T-nj+1);
Disc          = [Discl(101:-1:2);Disc];

%%
Ucs=[ ];  logCjs=[ ];

for i=1:ne
  s           = int2str(i);
  c           = cc(1:end,(i-1)*nj+1:i*nj);
  l           = lc(1:end,(i-1)*nj+1:i*nj);
  for j=2:nj;
    Cj          = diag(c(1:nj-j+1,j:nj));        
    Lj          = diag(l(1:nj-j+1,j:nj));        
    Ujlagc(j,1) = sum((beta.^[j-1:nj-1]').* ...
                  cumprod([diag(ps(1:nj-j+1,j-1:nj-1))]).* ...
                  (log(Cj)+ gam*(1-Lj).^(1-zeta(i))./(1-zeta(i))   ));
  
   loglagCj(j,1)= sum(Cj);    % sum(log(Cj));  
   loglagCj(j,1)= sum((beta.^[j-1:nj-1]').* ...
                               cumprod([diag(ps(1:nj-j+1,j-1:nj-1))]).* ...
                                   log(Cj) );   
  
  end;
  for t=1:T-nj+1;
    Cj          = diag(c(t:t+nj-1,:));        
    Lj          = diag(l(t:t+nj-1,:));        
    Ujc(t,1)    = sum((beta.^[0:nj-1]').* ...
                  cumprod([1;diag(ps(t:t+nj-1,1:nj-1))]).* ...
                  (log(Cj)+gam*(1-Lj).^(1-zeta(i))./(1-zeta(i))  ));
     logCj(t,1)= sum(Cj);     %sum(log(Cj));  
     logCj(t,1)= sum((beta.^[0:nj-1]').* ...
                        cumprod([1;diag(ps(t:t+nj-1,1:nj-1))]).* ...
                            log(Cj));   
  end;
  Ujc(T-nj+2:T)  = ones(nj-1,1)*Ujc(T-nj+1);
 
  eval(['Ujc',s,'=Ujc;']);  
  eval(['Uc',s,'= [Ujlagc(101:-1:2);Ujc];'])
   eval(['Ucs=[Ucs,Uc',s,'];']) %% add 2021/07/25
   
    logCj(T-nj+2:T)  = ones(nj-1,1)*logCj(T-nj+1);
    eval(['logCj',s,'=logCj;']);
    eval(['logC',s,'= [ loglagCj(101:-1:2) ; logCj ];'])
     eval(['logCjs=[logCjs,logC',s,'];'])
end;


%% New Policy 
 T = T+1;

cf            = zeros(nh,T);
lf            = zeros(nh,T);
ntx           = newpolt(23);
ntr           = newpolt(24);
% i1            = 503318;
% i1 = 504759+960+1;
% 
% cf(:)         = newpolt(i1+1:i1+nh*T);
% lf(:)         = newpolt(i1+1+nh*T:i1+2*nh*T);

i         = 26+nh;
i         = i+T;
i         = i+4*T*nh;
i         = i+2*nh;
i         = i+9*T;
i         = i+3*T*(ntx+ntr);
i         = i+23*T;
 cf(:)      =  newpolt(i+1+T*nh:i+2*T*nh);
 lf(:)      =  newpolt(i+1+2*T*nh:i+3*T*nh);

 %%
 cf(cf==0)=0.01;
 %%

cf            = cf';
lf            = lf';

cf = cf(2:end,:);
lf = lf(2:end,:);

if  retire_fg==1
lf(:,46:50)=0.001*ones(239,5);
lf(:,46+101:50+101)=0.001*ones(239,5);
lf(:,46+2*101:50+2*101)=0.001*ones(239,5);
lf(:,46+3*101:50+3*101)=0.001*ones(239,5);
end

T = T-1;

%%å¯ópä÷êîèCê≥ïKóv
Ufs = []; logCfs=[ ];
for i=1:ne
  s           = int2str(i);
  c           = cf(:,(i-1)*nj+1:i*nj);
  l           = lf(:,(i-1)*nj+1:i*nj);
  for j=2:nj;
    Cj          = diag(c(1:nj-j+1,j:nj));        
    Lj          = diag(l(1:nj-j+1,j:nj));   
    
    %  Moving average
%     type = 'linear';
%     windowSize = 3; % à⁄ìÆïΩãœÇÃä˙ä‘
%     Cj = movavg(Cj,type,windowSize);  % Moving Avara
%     Lj = movavg(Lj,type,windowSize);
    % ========================================
    
    Ujlagf(j,1) = sum((beta.^[j-1:nj-1]').* ...
                  cumprod([diag(ps(1:nj-j+1,j-1:nj-1))]).* ...
                  (log(Cj)+gam*(1-Lj).^(1-zeta(i))./(1-zeta(i))  ));
     loglagCj(j,1)= sum(Cj);    % sum(log(Cj));    
     loglagCj(j,1)= sum((beta.^[j-1:nj-1]').* ...
                               cumprod([diag(ps(1:nj-j+1,j-1:nj-1))]).* ...
                                     log(Cj) );   
              
  end;
  for t=1:T-nj+1;
    Cj          = diag(c(t:t+nj-1,:));        
    Lj          = diag(l(t:t+nj-1,:));    
    
    %  Moving average
%     type = 'linear';
%     windowSize = 5; % à⁄ìÆïΩãœÇÃä˙ä‘
%     Cj = movavg(Cj,type,windowSize);  % Moving Avara
%     Lj = movavg(Lj,type,windowSize);
    % ========================================
    
    Ujf(t,1)    = sum((beta.^[0:nj-1]').* ...
                  cumprod([1;diag(ps(t:t+nj-1,1:nj-1))]).* ...
                  (log(Cj)+gam*(1-Lj).^(1-zeta(i))./(1-zeta(i))  ));
       logCj(t,1)= sum(Cj);     %sum(log(Cj));     
       
       logCj(t,1)= sum((beta.^[0:nj-1]').* ...
                                    cumprod([1;diag(ps(t:t+nj-1,1:nj-1))]).* ...
                                  log(Cj));   
  end;
  Ujf(T-nj+2:T)  = ones(nj-1,1)*Ujf(T-nj+1);
  eval(['Ujf',s,'=Ujf;']);
  eval(['Uf',s,'= [Ujlagf(101:-1:2);Ujf];'])
   eval(['Ufs=[Ufs,Uf',s,'];']) %% add 2021/07/25
   
     logCj(T-nj+2:T)  = ones(nj-1,1)*logCj(T-nj+1);
    eval(['logCj',s,'=logCj;']);
    eval(['logC',s,'= [ loglagCj(101:-1:2) ; logCj ];'])
     eval(['logCfs=[logCfs,logC',s,'];'])
end;

% a=Ujlagf(101:-1:2);

% save('uf.mat','Ujf1','Uf1','Ujc1','Uc1','c','l')

WG = [];  Util=[];

for i=1:ne;
  s           = int2str(i);
  eval(['WG',s,'= (exp((Uf',s,'-Uc',s,')./Disc)-1)*100;'])  
    %Uf = new policy's Utility, Uc= current policy's Utility
  eval(['WG=[WG,WG',s,'];'])
  
%    eval(['Util_1',s,'= (exp((Uc',s,'-Uc1)./Disc))*100;'])   % current  policy
   eval(['Util_1',s,'= (exp((Uf',s,'-Uf1)./Disc))*100;'])   %  new policy
 %   eval(['Util_1',s,'= (Uc',s,'-Uc',s, ')*0 ;' ])  
    %Uf = new policy's Utility, Uc= current policy's Utility
  eval(['Util=[Util, Util_1',s,'];'])
  
end;


