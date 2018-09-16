%Plot um die Konvergenz bei feinerer Diskretisierung zu visualisieren:
figure1 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],real(DisplacementPorousNode1(:,1)),'xk')
plot([1:1:20],real(DisplacementPorousNode1(:,2)),'xb')
plot([1:1:20],real(DisplacementPorousNode1(:,3)),'.r')
xlabel("Lx/dx")
ylabel("Re(u_s_y(A))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

figure2 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],real(DisplacementPorousNode2(:,1)),'xk')
plot([1:1:20],real(DisplacementPorousNode2(:,2)),'xb')
plot([1:1:20],real(DisplacementPorousNode2(:,3)),'.r')
xlabel("Lx/dx")
ylabel("Re(u_s_y(B))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

figure3 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],real(DisplacementPorousNode3(:,1)),'xk')
plot([1:1:20],real(DisplacementPorousNode3(:,2)),'xb')
plot([1:1:20],real(DisplacementPorousNode3(:,3)),'.r')
xlabel("Lx/dx")
ylabel("Re(u_s_y(C))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

figure4 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],real(DisplacementPorousNode4(:,1)),'xk')
plot([1:1:20],real(DisplacementPorousNode4(:,2)),'xb')
plot([1:1:20],real(DisplacementPorousNode4(:,3)),'.r')
xlabel("Lx/dx")
ylabel("Re(u_s_y(D))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")


% Plot der Rechnenzeiten der Elemente abhängig von der Feinheit der Diskretisierung
figure4 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],real(DisplacementPorousTime(:,1)),'xk')
plot([1:1:20],real(DisplacementPorousTime(:,2)),'xb')
plot([1:1:20],real(DisplacementPorousTime(:,3)),'.r')
xlabel("Lx/dx")
ylabel("t in [s]")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")


