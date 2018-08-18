clc 
clear all
close all

DENSITY_S = 30;
DENSITY_F = 1.21;
POROSITY = 0.0001;
OMEGA = 100;
NUMBER_GAUSS_POINT = 2;

u=0;
for i = POROSITY:-0.00001:0.00001
u=u+1;
Differences(u,1)=i;
Differences(u,[2 3 4]) = ModelTest(i, DENSITY_S, DENSITY_F, OMEGA, NUMBER_GAUSS_POINT);
end