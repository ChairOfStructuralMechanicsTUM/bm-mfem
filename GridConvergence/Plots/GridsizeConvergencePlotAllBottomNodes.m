
figure1 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],real(DisplacementAllardSolid([1:1:91],1)),'*k')
plot([1:1:91],real(DisplacementAtallaSolid([1:1:91],1)),'xb')
plot([1:1:91],real(DisplacementDazelSolid([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("Re(u_s_y)")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")

figure2 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],imag(DisplacementAllardSolid([1:1:91],1)),'*k')
% plot([1:1:91],imag(DisplacementAtallaSolid([1:1:91],1)),'xb')
% plot([1:1:91],imag(DisplacementDazelSolid([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("Im(u_s_y)")
legend("(u_s,u_f)-Formulierung")

figure1 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:91],angle(DisplacementAllardSolid([1:1:91],1)),'*k')
plot([1:1:91],angle(DisplacementAtallaSolid([1:1:91],1)),'xb')
plot([1:1:91],angle(DisplacementDazelSolid([1:1:91],1)),'.r')
xlabel("Knoten-Id")
ylabel("\phi(u_s_y)")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")