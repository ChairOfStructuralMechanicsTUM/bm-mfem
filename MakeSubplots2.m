figure(1);     %<- Erste Figure die zu kopieren ist aktivieren
ax1 = gca;   %<- aktiven axishandle zuweisen
figure(2);
ax2 = gca;
figure(3);
ax3 = gca;

figure(4)       %<-Neue Figure
ax4 = subplot(1,3,1);      %<- Handle des Subplots
copyobj(allchild(ax1), ax4);

xlabel('Phase Im(k)')
ylabel('Frequenz f')
xlim([0 pi])
ylim([0 1.5e4]) 

ax5 = subplot(1,3,2);
copyobj(allchild(ax2), ax5); 

xlabel('Phase Im(k)')
ylabel('Frequenz f')
xlim([0 pi])
ylim([0 1.5e4]) 

ax6 = subplot(1,3,3);
copyobj(allchild(ax3), ax6); 

xlabel('Phase Im(k)')
ylabel('Frequenz f')
xlim([0 pi])
ylim([0 1.5e4]) 