
%Plot um die Konvergenz bei feinerer Diskretisierung zu visualisieren:

figure1 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,2),'xk')
%plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,4),'xk')
%plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,6),'xk')
xlabel("Lx/dx")
ylabel("u_s_y am Punkt nx")

figure2 = figure("Name","Solid y-Displacement Differences between Formulations");
hold on
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,4)-...
    DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,2),'xr')
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,6)-...
    DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,2),'xb')
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,6)-...
    DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,4),'xk')
xlabel("Lx/dx")
ylabel("u_s_y(i)-u_s_y(j) am Punkt nx")
legend("Atalla-Allard","Dazel-Allard","Dazel-Atalla")

figure3 = figure("Name","Solid y-Displacements low frequency range");
hold on
plot([1:1:20],DisplacementPorousEndNodeOmega12_5_F10_Betrag_Real(:,2),'xk')
plot([1:1:20],DisplacementPorousEndNodeOmega25_F10_Betrag_Real(:,2),'xr')
plot([1:1:5],DisplacementPorousEndNodeOmega50_F10_Betrag_Real(:,2),'xb')
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,2),'.k')
plot([1:1:20],DisplacementPorousEndNodeOmega150_F10_Betrag_Real(:,2),'.r')
xlabel("Lx/dx")
ylabel("u_s_y am Punkt nx")
legend("Omega = 12.5 Hz","Omega = 25 Hz","Omega = 50 Hz","Omega = 100 Hz","Omega = 150 Hz")


figure4 = figure("Name","Solid y-Displacements range 12.5-1600 Hz");
hold on
%plot([1:1:20],DisplacementPorousEndNodeOmega12_5_F10_Betrag_Real(:,2),'xk')
%plot([1:1:20],DisplacementPorousEndNodeOmega25_F10_Betrag_Real(:,2),'xr')
%plot([1:1:5],DisplacementPorousEndNodeOmega50_F10_Betrag_Real(:,2),'xb')
plot([1:1:20],DisplacementPorousEndNodeOmega100_F10_Betrag_Real(:,2),'.k')
%plot([1:1:20],DisplacementPorousEndNodeOmega150_F10_Betrag_Real(:,2),'.k')
plot([1:1:20],DisplacementPorousEndNodeOmega200_F10_Betrag_Real(:,2),'.c')
plot([1:1:20],DisplacementPorousEndNodeOmega400_F10_Betrag_Real(:,2),'.b')
plot([1:1:20],DisplacementPorousEndNodeOmega800_F10_Betrag_Real(:,2),'.g')
plot([1:1:20],DisplacementPorousEndNodeOmega1600_F10_Betrag_Real(:,2),'.r')
xlabel("Lx/dx")
ylabel("u_s_y am Punkt nx")
legend("Omega = 100 Hz","Omega = 200 Hz","Omega = 400 Hz","Omega = 800 Hz","Omega = 1600 Hz")