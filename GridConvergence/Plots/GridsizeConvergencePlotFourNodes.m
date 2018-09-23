%Plot um die Konvergenz bei feinerer Diskretisierung zu visualisieren:

%Array mit der Gesamtanzahl an Elementen als Wert für die Netzfeinheit:
for i  = 1:20
    a(i)=i*i*5
end
figure1 = figure("Name","Solid y-Displacements");
hold on
plot(a(:),real(DisplacementPorousNode1(:,2)),'*k')
plot(a(:),real(DisplacementPorousNode1(:,4)),'xb')
plot(a(:),real(DisplacementPorousNode1(:,6)),'.r')
xlabel("n_g_e_s")
ylabel("Re(u_s_y(A))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")
<<<<<<< HEAD
ylim([-6e-3 0])
=======
ylim([-7*10E-4 0])

>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076

figure2 = figure("Name","Solid y-Displacements");
hold on
plot(a(:),real(DisplacementPorousNode2(:,2)),'*k')
plot(a(:),real(DisplacementPorousNode2(:,4)),'xb')
plot(a(:),real(DisplacementPorousNode2(:,6)),'.r')
xlabel("n_g_e_s")
ylabel("Re(u_s_y(B))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")
<<<<<<< HEAD
ylim([-6e-3 0])
=======
ylim([-7*10E-4 0])
>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076

figure3 = figure("Name","Solid y-Displacements");
hold on
plot(a(:),real(DisplacementPorousNode3(:,2)),'*k')
plot(a(:),real(DisplacementPorousNode3(:,4)),'xb')
plot(a(:),real(DisplacementPorousNode3(:,6)),'.r')
xlabel("n_g_e_s")
ylabel("Re(u_s_y(C))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")
<<<<<<< HEAD
ylim([-6e-3 0])
=======
ylim([-7*10E-4 0])
>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076

figure4 = figure("Name","Solid y-Displacements");
hold on
plot(a(:),real(DisplacementPorousNode4(:,2)),'*k')
plot(a(:),real(DisplacementPorousNode4(:,4)),'xb')
plot(a(:),real(DisplacementPorousNode4(:,6)),'.r')
xlabel("n_g_e_s")
ylabel("Re(u_s_y(D))")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")
<<<<<<< HEAD
ylim([-6e-3 0])
=======
ylim([-7*10E-4 0])
>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076

% Plot der Rechnenzeiten der Elemente abhängig von der Feinheit der Diskretisierung
figure4 = figure("Name","Solid y-Displacements");
hold on
plot(a(:),real(DisplacementPorousTime(:,1)),'*k')
plot(a(:),real(DisplacementPorousTime(:,2)),'xb')
plot(a(:),real(DisplacementPorousTime(:,3)),'.r')
xlabel("n_g_e_s")
ylabel("t in [s]")
legend("(u_s,u_f)-Formulierung","(u_s,p)-Formulierung","(u_s,u_t)-Formulierung")


