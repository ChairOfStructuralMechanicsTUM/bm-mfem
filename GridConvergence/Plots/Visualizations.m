%SOLID

%Auslesen von Knotenkoordinaten

xx = modelDazel.getAllNodes.getX();

yy = modelDazel.getAllNodes.getY();

%Auslesen von Knotenverschiebungen

ux = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uy = modelDazel.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

%Skalierung zur besseren Visualisierung der Ergebnisse (damit man was

%sieht)

scaling = 1;

%Berechnen der Phase bzgl. ux oder uy; 41 und 21 sind hier die Anzahl der

%Knoten in x- bzw. y-Richtung.

<<<<<<< HEAD
z = reshape(imag(DisplacementDazelSolid(:,4)),91,19);
=======
z = reshape(sqrt(real(ux).^2+real(uy).^2),91,19);
>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076

%z = reshape(angle(uy),41,21);

%Berechnen der Knotenkoordinaten im Verformten System

xxx = reshape(xx+scaling*real(ux.'),91,19);

yyy = reshape(yy+scaling*real(uy.'),91,19);

%Abbilden der Ergebnisse

figure()


subplot(2,1,1)

surf(xxx,yyy,z,'FaceColor','interp')

xlabel("x [m]")
ylabel("y [m]")
<<<<<<< HEAD
c = colorbar
%c.Limits = [-1e-5 4e-6] 
=======
colorbar
%c.Limits = [-pi pi] 
>>>>>>> 59a973af27af5aacae8d1aeee79d4d59927c9076
view(0,90)