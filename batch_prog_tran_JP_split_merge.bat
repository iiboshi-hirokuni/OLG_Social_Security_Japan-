@echo off


set YYYYMMDD=%DATE:/=%



rem cd C:\Users\iibos\Dropbox\fortran\ozaki\code_adj

set PATH=C:\Program Files (x86)\mingw-w64\i686-8.1.0-posix-dwarf-rt_v6-rev0\mingw32\bin

rem set PATH=C:\cygwin\bin



rem ------------------------------------------------------------------------------------------
echo.
echo.

set /p BUNKI_01= " split file ? (y/n) -> " 

if %BUNKI_01% == y ( goto MODU_01
) else ( goto STEP_0
)


:MODU_01


gfortran -o split_file   split_file.f90


split_file

del split_file.o
del split_file.exe


:STEP_0
rem ------------------------------------------------------------------------------------------
echo.
echo.

set /p BUNKI_02= " merge file ? (y/n) -> " 

if %BUNKI_02% == y ( goto MODU_02
) else ( goto STEP_1
)


:MODU_02


gfortran -o merge_file   merge_file.f90


merge_file

del merge_file.o
del merge_file.exe
