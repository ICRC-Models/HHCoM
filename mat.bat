SET MATLABROOT="C:\Program Files\MATLAB\R2017b"
PATH=%MATLABROOT%;%PATH%
START matlab.exe -r %1 -logfile H:\temp\logfile -nodesktop -minimize -nosplash
PAUSE