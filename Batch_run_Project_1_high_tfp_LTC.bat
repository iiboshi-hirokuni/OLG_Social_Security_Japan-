@echo off


set YYYYMMDD=%DATE:/=%



rem cd C:\Users\iibos\Dropbox\fortran\ozaki\code_adj

set PATH=C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin

rem set PATH=C:\cygwin\bin



rem ------------------------------------------------------------------------------------------
echo.
echo.


rem " Start  transition Calculation"


:STEP_1

rem ------------------------------------------------------------------------------------------
echo.
echo.

set /p File_Name= "Input File Name   -> " 

copy .\file_in\trantf_JP_%File_Name%.inp  .\file_in\trantf_XX.inp 



set /p BUNKI_1= "complie ? (y/n) -> " 

if %BUNKI_1% == y ( goto MODU_1
) else ( goto STEP_4
)


:MODU_1

 gfortran -c .\fortran_codes\trans_JP_hs_rate_high_tfp.f90     >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1

:STEP_2

rem ------------------------------------------------------------------------------------------
echo.
echo.

set /p BUNKI_2= "LINK ? (y/n) -> "  

if %BUNKI_2% == y ( goto MODU_2
) else ( goto STEP_4
)

:MODU_2


gfortran -o trans_JP_hs_rate_high_tfp trans_JP_hs_rate_high_tfp.o .\mpi\mpi.o .\mpi\libmsmpi.a   >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1

:STEP_3

rem ------------------------------------------------------------------------------------------
echo.
echo.

set /p BUNKI_3= "parallel ? (y/n) -> "  

if %BUNKI_3% == y ( goto MODU_3
) else ( goto STEP_4
)

:MODU_3

rem ------------------------------------------------------------------------------------------
echo.
echo.


set /p HEIRETU_NUM=" Input Number of parallel processing ! (from 12 to 64)  -> "



rem set /p tpf_rate= "Input TPF Growth Rate (%) X 10  -> " 


 set tfp_rate=36
 
 
copy .\tfp_growth_rate3_LTC.txt  .\tfp_growth_rate.txt  >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1

"C:\Program Files\Microsoft MPI\Bin\mpiexec.exe" -n %HEIRETU_NUM% trans_JP_hs_rate_high_tfp.exe   >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1

@echo on


rem tfp_rate=6


copy  .\file_out\trantf_JP.dat .\file_out\trantf_%File_Name%_grate%tfp_rate%_LTC.dat >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1

copy .\file_out\trantf_JP.nxt .\file_out\trantf_%File_Name%.nxt >> result_%File_Name%_high_TFP_%YYYYMMDD%.txt 2>&1




:STEP_4

del trans_JP_hs_rate_high_tfp.o
del trans_JP_hs_rate_high_tfp.exe
