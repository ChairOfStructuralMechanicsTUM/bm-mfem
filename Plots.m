
%Plot um die Konvergenz bei feinerer Diskretisierung zu visualisieren:
figure1 = figure("Name","Solid x-Displacements");
hold on
plot([1:1:size(DisplacementPorousEndNode(:,1),1)],DisplacementPorousEndNode(:,1),'xk')
plot([1:1:size(DisplacementPorousEndNode(:,1),1)],DisplacementPorousEndNode(:,3),'xb')
plot([1:1:size(DisplacementPorousEndNode(:,1),1)],DisplacementPorousEndNode(:,5),'xr')
xlabel("Diskretisierungsfaktor i")
ylabel("Verschiebung am Punkt nx")

figure2 = figure("Name","Solid y-Displacements");
hold on
plot([1:1:size(DisplacementPorousEndNode(:,1),1)],DisplacementPorousEndNode(:,2),'xk')
plot([1:1:size(DisplacementPorousEndNode(:,2),1)],DisplacementPorousEndNode(:,4),'xb')
plot([1:1:size(DisplacementPorousEndNode(:,2),1)],DisplacementPorousEndNode(:,6),'xr')
xlabel("Diskretisierungsfaktor i")
ylabel("Verschiebung am Punkt nx")

