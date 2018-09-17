figure1 = figure("Name","");
hold on
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,1),'*k')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,2),'xb')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,3),'.r')
plot([0.04:0.1:0.94],DisplacementPorousNode4_96_006(:,4),'-k')
xlabel("1- phi")
ylabel("Re(u_s_y(nx+1))")
legend("(u_s,u_f) - Formulierung","(u_s,p) - Formulierung","(u_s,u_t) - Formulierung","Homogenes Element")

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
