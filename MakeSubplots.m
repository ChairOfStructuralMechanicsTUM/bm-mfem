figure(1);     %<- Erste Figure die zu kopieren ist aktivieren
ax1 = gca;   %<- aktiven axishandle zuweisen
figure(2);
ax2 = gca;
figure(3)       %<-Neue Figure
ax3 = subplot(1,2,1);      %<- Handle des Subplots
copyobj(allchild(ax1), ax3);


ax4 = subplot(1,2,2);
copyobj(allchild(ax2), ax4); 

xlabel('Phase k')
ylabel('Frequenz f')
xlim([0 pi])
ylim([0 1.5e4]) 