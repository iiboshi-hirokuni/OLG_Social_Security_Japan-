# tran_JP_hs_rate_v7の実行方法

[toc]



## Step 1:  tfp_growth_rate3.txtの値の変更

- #### LTCの場合

　　　**tfp_growth_rate3_LTC.txt**の書き換え
　　　　+0.10    over age 75: Rate of increase of copay of **long term care**
　　　　+0.10    age 66 thru 75: copay of long term care

- #### MEDの場合

　　**tfp_growth_rate3_LTC.txt**の書き換え

​          0.10    over age 75: Rate of increase of Co pay of medical expenditure 
​         0.10    age 66 thru 75: Co pay of medical expenditure

- ### tfp_growth_rate3.txtの値の変更方法

##### 医療費　の自己負担率が10%上昇する場合　

0.003  tfpの成長率
0.10    over age 75: Rate of increase of Co pay of medical expenditure
0.10    age 66 thru 75: Co pay of medical expenditure
0.00    over age 75: Rate of increase of copay of long term care
0.00    age 66 thru 75: copay of long term care



## Step 2: インプットファイルの選択

#### RR50%の場合

trantf_JP_muprop_sec2_tfp1_gz07_gn_rr50.inp

#### それ以外

trantf_JP_muprop_sec2_tfp1_gz07_gn.inp



## Step 3: バッチファイルの選択

#### シナリオ1　現状維持

Batch_run_Project_1.bat

#### シナリオ2　2030年に定年延長

Batch_run_Project_1_retire_ext.bat

#### シナリオ3　2030年に　tfp成長率が0.6%に上昇

Batch_run_Project_1_high_tfp.bat




### バッチの入力イメージ



![batch](C:\Users\iiboshi\Downloads\batch.png)











