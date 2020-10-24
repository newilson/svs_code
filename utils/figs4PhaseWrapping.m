t = [0,1,2];
phi = [20,120,220];
phiwrap = [20,120,-140];

figure, plot(t,phi,'o',t,phiwrap,'*')
legend({'true','measured'})
ylim([-360 360]),xlim([-.1 2.1])
xlabel('time'), ylabel('phase')
print(gcf,'Fig1','-dpng','-r300');