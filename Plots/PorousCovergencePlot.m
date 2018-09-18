figure1 = figure("Name","");
hold on
<<<<<<< HEAD
plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,2),'*k')
% plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,4),'xb')
% plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,6),'.r')
plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,8),'-k')
legend("(u_s,u_f) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-12;1e-12]")
=======
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,1),'*k')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,2),'xb')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,3),'.r')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,u_f) - Formulierung","(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
>>>>>>> cbb824cb0bf6330a0bc383d69cf1349c690b7122

figure1 = figure("Name","");
hold on
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,1),'*k')
% plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,2),'xb')
% plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,3),'.r')
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,u_f) - Formulierung","Homogenes Element")


figure1 = figure("Name","");
hold on
<<<<<<< HEAD
% plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,2),'*k')
plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,4),'xb')
plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,6),'.r')
plot([0.999999999991:0.000000000001:0.999999999999],DisplacementPorousNonPorousEndNode_0000000000009_0000000000001(:,8),'-k')
legend("(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-12;1e-12]")


figure3 = figure("Name","");
hold on
plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,2),'*k')
% plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,4),'xb')
% plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,6),'.r')
plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,8),'-k')
legend("(u_s,u_f) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-8;1e-8]")

figure4 = figure("Name","");
hold on
% plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,2),'*k')
plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,4),'xb')
plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,6),'.r')
plot([0.99999991:0.00000001:0.99999999],DisplacementPorousNonPorousEndNode_000000009_000000001(:,8),'-k')
legend("(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-8;1e-8]")

figure5 = figure("Name","");
hold on
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,2),'*k')
% plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,4),'xb')
% plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,6),'.r')
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,8),'-k')
legend("(u_s,u_f) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-3;1e-3]")

figure6 = figure("Name","");
hold on
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,2),'*k')
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,4),'xb')
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,6),'.r')
plot([0.991:0.001:0.999],DisplacementPorousNonPorousEndNode_009_001(:,8),'-k')
legend("(u_s,u_f) - Formulierung","(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [9e-3;1e-3]")

figure7 = figure("Name","");
hold on
plot([0.04:0.1:0.94],DisplacementPorousNonPorousEndNode_96_06(:,2),'*k')
plot([0.04:0.1:0.94],DisplacementPorousNonPorousEndNode_96_06(:,4),'xb')
plot([0.04:0.1:0.94],DisplacementPorousNonPorousEndNode_96_06(:,6),'.r')
plot([0.04:0.1:0.94],DisplacementPorousNonPorousEndNode_96_06(:,8),'-k')
legend("(u_s,u_f) - Formulierung","(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
xlabel("1-\phi")
ylabel("Re(u_s_y(91))")
title("\phi = [0.96;0.06]")
=======
% plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,1),'*k')
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,2),'xb')
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,3),'.r')
plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")

figure1 = figure("Name","");
hold on
% plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,1),'*k')
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,2),'xb')
plot([0.991:0.001:0.999],DisplacementPorousNode4_009_001(:,3),'.r')
% plot([0.991:0.001:0.999],DisplacementPorousNode1_009_001(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,p) - Formulierung","(u_s,u_t) - Formulierung")


figure1 = figure("Name","");
hold on
 %plot([(1-9*10^(-10)):10^(-10):(1-1*10^(-10))],DisplacementPorousNode1_9e10_1e10(:,1),'*k')
 plot([(1-9*10^(-10)):10^(-10):(1-1*10^(-10))],DisplacementPorousNode1_9e10_1e10(:,2),'xb')
 plot([(1-9*10^(-10)):10^(-10):(1-1*10^(-10))],DisplacementPorousNode1_9e10_1e10(:,3),'.r')
 plot([(1-9*10^(-10)):10^(-10):(1-1*10^(-10))],DisplacementPorousNode1_9e10_1e10(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")
>>>>>>> cbb824cb0bf6330a0bc383d69cf1349c690b7122
