@echo off

set PYTHONFOLDER=C:\Python27

set IMARISFOLDER=C:\Program Files\Bitplane\Imaris x64 9.5.0
set PATH=%IMARISFOLDER%;%PATH%

set PYTHONPATH=%IMARISFOLDER%\XT\python2;%IMARISFOLDER%\XT\python2\private;%PYTHONPATH%

%PYTHONFOLDER%\python.exe %PYTHONFOLDER%\Lib\idlelib\idle.pyw