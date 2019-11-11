  M = dlmread('osilator_fundamental.dat');
  wl = M(:,1);
  Int = M(:,2)./max(M(:,2));
  plot(wl,Int)
  title('Spectrum')
  xlim([500 1200])
  ylim([0 1.2])
  yticks(0:1:1)
  xticks([600 700 800 900 1030 1200])
 
  
