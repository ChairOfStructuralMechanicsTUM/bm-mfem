% testAnsysInput

clear all
close all
clc

file='F:\bmFEM\bm-mfem\examples\3DBeam.txt';
ansysExecutable='C:\Program Files\ANSYS Inc\v171\ansys\bin\winx64\ANSYS171.exe';

ansysInput = AnsysInput(file);

data=ansysInput.readModel(ansysExecutable);

