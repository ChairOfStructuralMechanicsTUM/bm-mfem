function [ path ] = getPaths_changeme(type)
%GETPATHS_CHANGEME This file contains the paths to required executables
%   Copy this file and change all paths according to your system. Then
%   change the filename to "getPaths.m" and the function name to "getPaths".
switch type
    case 'ansys'
        path = 'C:/Program Files/ANSYS Inc/v171/ansys/bin/winx64/ANSYS171.exe';
    otherwise
        msg = ['GetPaths: A path with name \"', ...
            type, '\" is not defined'];
        e = MException('MATLAB:bm_mfem:undefinedPathName',msg);
        throw(e);
end
end

