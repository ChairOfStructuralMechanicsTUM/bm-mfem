%Laufzeitplots 
x = [100, 225, 400, 625, 900, 1225];
y1=[0.4946 1.1013 2.0366 3.4134 5.1736 7.4519];
plot(x,y1)
hold on 
y2=[0.6736 1.7356 4.5660 11.1733 28.0187 64.0644];
plot(x,y2)
legend('Laufzeit-FEM','Laufzeit-FETI-DP','Location','northwest')
xlabel('Knotenzahlen der Systeme')
ylabel('Laufzeit [s]')