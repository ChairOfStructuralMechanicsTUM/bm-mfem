function minLambda = WavelengthCheck(omega,alpha_inf,eta_f,rho_f,rho_s,Xi,Lambda_v,phi,gamma,P_0,...
    Pr,Lambda_t,eta_s,lambda,mu,Lx,Ly,nx,ny)

% G_J(OMEGA) = flow resistivity of air particles in the pores:
G_J = (1+(4*1i*omega*alpha_inf^2*eta_f*rho_f)/...
    (Xi^2*Lambda_v^2*phi^2))^(0.5);
% Viscous Drag accounts for viscous body forces interacting between solid and fluid phases,
% proportional to the relative velocity:
bF = Xi*phi^2*G_J;
% Calculate Bulk-Modulus K_f:
K_f = (gamma*P_0)/(gamma-(gamma-1)*...
    (1+(8*eta_f)/(1i*omega*Pr*Lambda_t^2*rho_f)*...
    (1+(1i*omega*Pr*Lambda_t^2*rho_f)/(16*eta_f))^(0.5))^(-1));
LameCoeff = (1+1i*eta_s)*lambda;
N = (1+1i*eta_s)*mu;
% Biot's Elasticity Coefficients:
R = phi*K_f;
%Q = (1-POROSITY)*K_f;
A = LameCoeff+(1-phi)^2/phi*K_f;
P = A+2*N;

rho_ss = (1-phi)* rho_s - 1i*bF/omega;
rho_ff = phi*P_0-1i*bF/omega;
rho_sf = 1i*bF/omega;

k1 = sqrt(omega^2*(rho_ss/P + rho_sf^2/(rho_ss*R-rho_ff*P)));
k2 = sqrt(omega^2*(rho_ff/R-rho_sf^2/(rho_ss*R-rho_ff*P)));
k3 = sqrt(omega^2/N*(1-phi)*rho_s);

lambda1=1/k1;
lambda2=1/k2;
lambda3=1/k3;

minLambda = min([lambda1,lambda2,lambda3]);
% if 10<min([lambda1,lambda2,lambda3])/max([Lx/nx,Ly/ny])
%     fprintf("error: Wavelength-Ratio: lambda < (max(dx,dy))/10 \n ")
%     return
% end


end