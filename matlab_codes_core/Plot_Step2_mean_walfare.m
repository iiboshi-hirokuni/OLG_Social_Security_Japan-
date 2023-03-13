% 
% 
%  Calculation of Walfare for Retierment, Current & Future Generations  
%  
% 
% 

%close all

StartAge = 20;   % From 100 Age 
EndAge  = 200;  % Upto -60 Age

WG_Table ={};   Mean_WG_S1={};   Mean_WG_S2={} ;  

WG_Table{1} = WG(1+StartAge:EndAge,:);  %  All Generations
WG_Table{2} = WG(1+StartAge:1+StartAge+(100-65),:);  % Retirement from 100 Age to 65 Age
WG_Table{3} = WG(1+StartAge+(100-64):1+StartAge+(100-20),:);  % Current From 64 Age to 20 Age
WG_Table{4} = WG(1+StartAge+(100-19):EndAge,:);  % Future from 19 Age to -60 age

varNames =   {  'Agents  '    , 'Total',  'Retire', 'Current','Future' };
Agent_name =  [ "male-regular";  "female-regular"; "male-contingent"; "female-contingent" ];        

if Scenario ==1
       disp( ' Make  Table of Senario 1   '   ) 
     
        WG_S1 = WG_Table;
        for  i = 1:4
          Mean_WG_S1{i} = mean(WG_S1{i});
       end   
        
       Tab=table(Agent_name ,Mean_WG_S1{1}', ...
                           Mean_WG_S1{2}', Mean_WG_S1{3}', Mean_WG_S1{4}',...
                           'VariableNames',varNames);
       
       writetable(Tab,'Table-mean-Senario1.csv')
        
else  %Scenario ==2
      disp( ' Make  Table of Senario 2   '   ) 
    
       WG_S2 = WG_Table;
       for  i = 1:4
          Mean_WG_S2{i} = mean(WG_S2{i});
       end   
        
       Tab=table(Agent_name ,Mean_WG_S2{1}', ...
                           Mean_WG_S2{2}', Mean_WG_S2{3}', Mean_WG_S2{4}',...
                           'VariableNames',varNames);
       
       writetable(Tab,'Table-mean-Senario2.csv')
       
end

