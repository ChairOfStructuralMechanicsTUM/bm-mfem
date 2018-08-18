clc 
clear all
close all

DENSITY_S = 30;
DENSITY_F = 1.21;
POROSITY = 0.96;
OMEGA = 100;
NUMBER_GAUSS_POINT = 2;

ModelTest(POROSITY, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT)
% u=1;
% Differences(u,1)=i;
% Differences(u,[2 3 4]) = ModelTest(i, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT);
