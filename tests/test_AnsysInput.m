% testAnsysInput

clear all
close all
clc

file='F:\ALMA_simulations\Tools\SVN_Dispersion\Models\3DBeam.txt';
ansysExecutable='C:\Program Files\ANSYS Inc\v171\ansys\bin\winx64\ANSYS171.exe';

ansysInput = AnsysInput(file);

data=ansysInput.readModel(ansysExecutable);

