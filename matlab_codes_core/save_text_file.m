
model_name = 'BASELINE';
% model_name = 'TFP';

est_date = datestr(date);         
result_name = ['./output/ESTIMATE_', char(model_name), '_' ,est_date, '.txt'];          

fileID = fopen(result_name,'w');

% fprintf(fileID

fprintf(fileID,'Results_of_OLG_Model_with_Intangible_Capital \n');
fprintf(fileID,'--------------------------------------------\n');
fprintf(fileID,'  Parameters \n');
fprintf(fileID,'  ---------- \n');
fprintf(fileID,sprintf('    Asset_grid_lower_bound               %8.3f \n',alb));
fprintf(fileID,sprintf('    Asset_grid_upper_bound               %8.3f \n',aub));
fprintf(fileID,sprintf('    Corporate_share_of_output            %8.3f \n',alpha));
fprintf(fileID,sprintf('    Discount_factor                      %8.3f \n',beta));
fprintf(fileID,sprintf('    Tangible_depreciation,_sector_1      %8.3f \n',delt1));
fprintf(fileID,sprintf('    Intangible_depreciation,_sector_1    %8.3f \n',deli1));
fprintf(fileID,sprintf('    Tangible_depreciation,_sector_2      %8.3f \n',delt2));
fprintf(fileID,sprintf('    Intangible_depreciation,_sector_2    %8.3f \n',deli2));
fprintf(fileID,sprintf('    Weight_on_leisure_in_utility         %8.3f \n',gam));
fprintf(fileID,sprintf('    Growth_of_technology                 %8.3f \n',gz*100));
fprintf(fileID,sprintf('    Tangible_capital_share,_sector_1     %8.3f \n',thett1));
fprintf(fileID,sprintf('    Intangible_capital_share,_sector_1   %8.3f \n',theti1));
fprintf(fileID,sprintf('    Tangible_capital_share,_sector_2     %8.3f \n',thett2));
fprintf(fileID,sprintf('    Intangible_capital_share,_sector_2   %8.3f \n',theti2));
fprintf(fileID,sprintf('    Age_of_first_job                     %8g \n',21));
fprintf(fileID,sprintf('    Age_when_retired                         \n'));
for i=1:ne
fprintf(fileID,sprintf(['      Group_',int2str(i), ...
              '                            %8g \n '],Jr(i)+20));
end
fprintf(fileID,sprintf('    Maximum_age                          %8g \n ',nj+20));
fprintf(fileID,sprintf('    Degree_of_annuitization              %8g \n',lam*100));
fprintf(fileID,'  \n ');
fprintf(fileID,'  \n ');

fprintf(fileID,'                                                   Period_Averages  \n');
fprintf(fileID,'  ----------------                --------------------------------------------------       \n   ') ;
fprintf(fileID,'  Exogenous_Inputs              "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,'  ----------------                ------------------------------------------------------------ \n');
fprintf(fileID,sprintf('    Government_borrowing       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  \n',AVG(:,1)));
fprintf(fileID,sprintf('    Government_spending        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  \n',AVG(:,2)));
fprintf(fileID,'    Tax_rates \n');
fprintf(fileID,sprintf('      Consumption              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,3)));
fprintf(fileID,sprintf('      Distributions            %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,4)));
fprintf(fileID,sprintf('      Labor                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,5)));
fprintf(fileID,sprintf('      Corporate_profits        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,6)));
fprintf(fileID,sprintf('      Noncorporate_profts      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,7)));
fprintf(fileID,sprintf('    TFP_parameters             %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,8)));
fprintf(fileID,sprintf('    Transfers_to_retirees      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,71+3*ne)));
for i=1:ne;
fprintf(fileID,sprintf(['      Group',int2str(i), ...
              '                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n '],AVG(:,71+3*ne+i)));
end
fprintf(fileID,sprintf('    Fraction_working           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n ',AVG(:,70+2*ne)));
for i=1:ne;
fprintf(fileID,sprintf(['      Group',int2str(i), ...
              '                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n '],AVG(:,70+2*ne+i)));
end
fprintf(fileID,' \n');
fprintf(fileID,'  Endogenous_Variables,_averaged   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,'  ------------------------------ \n');
fprintf(fileID,sprintf('    Interest_rate              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  \n',AVG(:,9)));
fprintf(fileID,sprintf('    Wage_rate                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  \n',AVG(:,10)));
fprintf(fileID,sprintf('    HH_transfer                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f  \n',AVG(:,11)));
fprintf(fileID,' \n');
fprintf(fileID,'  National_Accounts_(rel._to_GDP) \n ');
fprintf(fileID,'  ------------------------------- \n ');
fprintf(fileID,'  *Gross_National_Product/Income   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,' \n')  ;
fprintf(fileID,sprintf('    Consumption_(C)            %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,12)));
fprintf(fileID,sprintf('    Tangible_investment_(XT)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,13)));
fprintf(fileID,sprintf('      Corporate_(XT1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,14)));
fprintf(fileID,sprintf('      Noncorporate_(XT1)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,15)));
fprintf(fileID,sprintf('    Government_Spending_(G)    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,16)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    GDP_(Y)                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',ones(1,7)));
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Depreciation_(delT_KT)     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,17)));
fprintf(fileID,sprintf('      Corporate_(delT_KT1)     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,18)));
fprintf(fileID,sprintf('      Noncorporate_(delT_KT2)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,19)));
fprintf(fileID,sprintf('    Labor_compensation_(wL)    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,20)));
        fprintf(fileID,'    Profits,_sum_pi_Yi-wLi-depTi_KTi-XIi \n ');
fprintf(fileID,sprintf('      Corporate                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,21)));
fprintf(fileID,sprintf('        Taxes                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,22)));
fprintf(fileID,sprintf('        Net_dividends          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,23)));
fprintf(fileID,sprintf('        Undistributed          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,24)));
fprintf(fileID,sprintf('      Noncorporate             %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,25)));
fprintf(fileID,sprintf('        Taxes                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,26)));
fprintf(fileID,sprintf('        Pre-tax_distributions  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,27)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    GDI_=_GDP                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',ones(1,7)));
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Addenda:_intangible_(XI)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,28)));

fprintf(fileID,' \n');
fprintf(fileID,'  *Government_Expenditures/Receipts  "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,' \n')  ;
fprintf(fileID,sprintf('    Consumption_of_G&S         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,29)));
fprintf(fileID,sprintf('    Transfers                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,30)));
fprintf(fileID,sprintf('      To_retired               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,31)));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
             '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,71+4*ne+i)));
end;
fprintf(fileID,sprintf('      To_non-retired           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,32)));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
             '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,71+5*ne+i)));
end;
fprintf(fileID,sprintf('    Interest_paid              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,34)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    Expenditures               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,35)));
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Tax_revenues               %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,36)));
fprintf(fileID,sprintf('      Consumption              %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,37)));
fprintf(fileID,sprintf('      Labor                    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,38)));
for i=1:ne;
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG2(:,i)));
end;
fprintf(fileID,sprintf('      Dividends                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,39)));
fprintf(fileID,sprintf('      C-Profits                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,40)));
fprintf(fileID,sprintf('      U-Profits                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,41)));
fprintf(fileID,sprintf('    Interest_received          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,42)));
fprintf(fileID,sprintf('    Accidental_bequests        %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,43)));
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Receipts                   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,44)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    Deficit=Expend.-Receipts   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,45)));
fprintf(fileID,' \n');
fprintf(fileID,'  Fixed_Assets_(rel._to_GDP)   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,'  -------------------------- \n ');
fprintf(fileID,sprintf('    Tangible_capital_(KT'')     %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,46)));
fprintf(fileID,sprintf('      Corporate_(KT1'')         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,47)));
fprintf(fileID,sprintf('      Noncorporate_(KT2'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,48)));
fprintf(fileID,sprintf('    Addenda:_intangible_(KI'')  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,49)));
fprintf(fileID,sprintf('      Corporate_(KI1'')         %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,50)));
fprintf(fileID,sprintf('      Noncorporate_(KI2'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,51)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    Total                      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,52)));
fprintf(fileID,' \n');
fprintf(fileID,'  Flow_of_Funds_(rel._to_GDP) \n ');
fprintf(fileID,'  --------------------------- \n ');
fprintf(fileID,'  *Household_Balance_Sheet   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Tangible_capital_(vT)      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,53)));
fprintf(fileID,sprintf('      Corporate_(vT1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,54)));
fprintf(fileID,sprintf('      Noncorporate_(vT2)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,55)));
fprintf(fileID,sprintf('    Government_Debt_(B'')       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,56)));
fprintf(fileID,sprintf('    Addenda:_intangible_(vI)   %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,57)));
fprintf(fileID,sprintf('      Corporate_(vI1)          %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,58)));
fprintf(fileID,sprintf('      Noncorporate_(vI2)       %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,59)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    Net_Worth                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,60)));
fprintf(fileID,' \n');
fprintf(fileID,'  *Government_Balance_Sheet   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,' \n');
fprintf(fileID,sprintf('    Tangible_capital_(KT''-vT)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,61)));
fprintf(fileID,sprintf('    Government_Debt_(-B'')      %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,62)));
fprintf(fileID,sprintf('    Addenda:_intang._(KI''-vI)  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,63)));
fprintf(fileID,     '                                  =====    =====    =====    =====    =====    =====    ===== \n');
fprintf(fileID,sprintf('    Net_Worth                  %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,64)));
fprintf(fileID,' \n');
fprintf(fileID,'  Welfare \n');
fprintf(fileID,'  ------- \n');
fprintf(fileID,'  *Individual_consumption   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,sprintf('      Weighted_sum_mu_j*c_j    %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,66)));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,71+6*ne+i)));
end;
fprintf(fileID,'  *Individual_labor_supply   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,sprintf('      Weighted_sum_mu_j*l_j*e_j%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,68)));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,71+7*ne+i)));
end;
fprintf(fileID,'  *Values_for_computing_welfare   "1-10"    "11-20"    21-30    31-60    61-89   90-200  BG_Path \n');
fprintf(fileID,sprintf('      Lifetime_utility_(U) \n'));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,69+i)));
end;
fprintf(fileID,sprintf('      Discount_(sum_beta^j_ps(t,j)) \n '));
for i=1:ne
fprintf(fileID,sprintf(['        Group_',int2str(i), ...
              '                %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n'],AVG(:,69+ne+i)));
end;
fprintf(fileID,' \n');
fprintf(fileID,'    % Welfare gain=(exp((U1-U2)/Discount)-1)*100 \n ');
fprintf(fileID,' \n');
fprintf(fileID,'  --------------------------- \n');
fprintf(fileID,' \n');
fprintf(fileID,'  NOTES: \n');
fprintf(fileID,sprintf('    (1)_level_of_GDP           %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f \n',AVG(:,69)));
fprintf(fileID,' \n');


