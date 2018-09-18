%SOLID

%Auslesen von Knotenkoordinaten

xx = modelAllard.getAllNodes.getX();

yy = modelAllard.getAllNodes.getY();

%Auslesen von Knotenverschiebungen

ux = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_X');

uy = modelAllard.getAllNodes.getDofValue('DISPLACEMENT_SOLID_Y');

%Skalierung zur besseren Visualisierung der Ergebnisse (damit man was

%sieht)

scaling = 1;

%Berechnen der Phase bzgl. ux oder uy; 41 und 21 sind hier die Anzahl der

%Knoten in x- bzw. y-Richtung.

z = reshape(imag(ux),91,19);

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
colorbar

view(0,90)