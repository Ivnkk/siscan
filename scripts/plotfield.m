

Np = length(f0);
% fs = 5;
% dt = 1/fs;
dt = mean((f0(end)- f0(1))./f0.^2);
% t =(-Np/2:Np/2-1)*dt/Np;
t = linspace(-dt*(Np-1)/2,dt*(Np-1)/2,Np);

Ret = ifft(field);
Ret = ifftshift(Ret);
phase = unwrap(angle(Ret));
Inten = abs(Ret).^2./max(abs(Ret).^2);

plot(t,Inten)
yyaxis left
title ('Temporal intensity profile')
ylabel('Intensity, a.u.')
xlabel('time, fs')
ylim([0 1.2])
yyaxis right
plot(t,phase)
ylabel('Temporal phase, rad')


t3 = (-1 + dt/Np):(dt/Np):(1-dt/Np);