function minLambda = WavelengthCheck(OMEGA,TORTUOSITY,ETA_F,DENSITY_F,DENSITY_S,STATIC_FLOW_RESISTIVITY,VISCOUS_CHARACT_LENGTH,POROSITY,HEAT_CAPACITY_RATIO,PRESSURE_0,...
    PRANDL_NUMBER,THERMAL_CHARACT_LENGTH,ETA_S,LAMBDA,MUE,Lx,Ly,nx,ny)

% G_J(OMEGA) = flow resistivity of air particles in the pores:
AirFlowResistivity = (1+(4*1i*OMEGA*TORTUOSITY^2*ETA_F*DENSITY_F)/...
    (STATIC_FLOW_RESISTIVITY^2*VISCOUS_CHARACT_LENGTH^2*POROSITY^2))^(0.5);

% Viscous Drag accounts for viscous body forces interacting between solid and fluid phases,
% proportional to the relative velocity:
bF = STATIC_FLOW_RESISTIVITY*POROSITY^2*AirFlowResistivity;

% Calculate Bulk-Modulus K_f:
K_f = (HEAT_CAPACITY_RATIO*PRESSURE_0)/(HEAT_CAPACITY_RATIO-(HEAT_CAPACITY_RATIO-1)*...
    (1+(8*ETA_F)/(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)*...
    (1+(1i*OMEGA*PRANDL_NUMBER*THERMAL_CHARACT_LENGTH^2*DENSITY_F)/(16*ETA_F))^(0.5))^(-1));

LameCoeff = (1+1i*ETA_S)*LAMBDA;
N = (1+1i*ETA_S)*MUE;

% Biot's Elasticity Coefficients:
R = POROSITY*K_f;
%Q = (1-POROSITY)*K_f;
A = LameCoeff+(1-POROSITY)^2/POROSITY*K_f;
P = A+2*N;

rho_ss = (1-POROSITY)* DENSITY_S - 1i*bF/OMEGA;
rho_ff = POROSITY*PRESSURE_0-1i*bF/OMEGA;
rho_sf = 1i*bF/OMEGA;

k1 = sqrt(OMEGA^2*(rho_ss/P + rho_sf^2/(rho_ss*R-rho_ff*P)));
k2 = sqrt(OMEGA^2*(rho_ff/R-rho_sf^2/(rho_ss*R-rho_ff*P)));
k3 = sqrt(OMEGA^2/N*(1-POROSITY)*DENSITY_S);

lambda1=1/k1;
lambda2=1/k2;
lambda3=1/k3;

minLambda = min([lambda1,lambda2,lambda3]);
% if 10<min([lambda1,lambda2,lambda3])/max([Lx/nx,Ly/ny])
%     fprintf("error: Wavelength-Ratio: lambda < (max(dx,dy))/10 \n ")
%     return
% end


end