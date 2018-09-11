% Unteren Rechts Differenzen
hold on
plot([nx-9:nx+1],real( DisplacementDifferencePorous([nx-9:nx+1],1)),'xk')
plot([nx-9:nx+1],real( DisplacementDifferencePorous([nx-9:nx+1],2)),'.b')
plot([nx-9:nx+1],real( DisplacementDifferencePorous([nx-9:nx+1],3)),'xr')
legend("Allard-Atalla","Allard-Dazel","Dazel-Atalla")
xlabel("Node-Id")
ylabel("u_s_y(Node-Id)")
%title("Displacement Differences Node by Node")


% Unteren Rechts Absoulut
hold on
plot([nx-9:nx+1],real( DisplacementPorous([nx-9:nx+1],1)),'xk')
plot([nx-9:nx+1],real( DisplacementPorous([nx-9:nx+1],2)),'.b')
plot([nx-9:nx+1],real( DisplacementPorous([nx-9:nx+1],3)),'xr')
legend("Allard-Atalla","Allard-Dazel","Dazel-Atalla")
xlabel("Node-Id")
ylabel("u_s_y(Node-Id)")
%title("Displacement Differences Node by Node")


% Unterer Rand Differenzen
hold on
plot([1:nx+1],real( DisplacementDifferencePorous([1:nx+1],1)),'xk')
plot([1:nx+1],real( DisplacementDifferencePorous([1:nx+1],2)),'.b')
plot([1:nx+1],real( DisplacementDifferencePorous([1:nx+1],3)),'xr')
legend("Allard-Atalla","Allard-Dazel","Dazel-Atalla")
xlabel("Node-Id")
ylabel("u_s_y(Node-Id)")
%title("Displacement Differences Node by Node")

% Unterer Rand Absolut
hold on
plot([1:nx+1],real( DisplacementPorous([1:nx+1],1)),'xk')
plot([1:nx+1],real( DisplacementPorous([1:nx+1],2)),'.b')
plot([1:nx+1],real( DisplacementPorous([1:nx+1],3)),'xr')
legend("Allard-Atalla","Allard-Dazel","Dazel-Atalla")
xlabel("Node-Id")
ylabel("u_s_y(Node-Id)")
%title("Displacement Differences Node by Node")

% Oberer Rand
hold on
plot([(ny-1)*(nx+1):ny*(nx+1)],real( DisplacementDifferencePorous([(ny-1)*(nx+1):ny*(nx+1)],1)),'xk')
plot([(ny-1)*(nx+1):ny*(nx+1)],real( DisplacementDifferencePorous([(ny-1)*(nx+1):ny*(nx+1)],2)),'.b')
plot([(ny-1)*(nx+1):ny*(nx+1)],real( DisplacementDifferencePorous([(ny-1)*(nx+1):ny*(nx+1)],3)),'xr')
legend("Allard-Atalla","Allard-Dazel","Dazel-Atalla")
xlabel("Node-Id")
ylabel("u_s_y(Node-Id)")
%title("Displacement Differences Node by Node")