function [Inten,phase,t] = getField_t(scan,wf)
%GETFIELD_T Summary of this function goes here
%   Detailed explanation goes here
Np = length(wf);
dt = mean(2*pi*(wf(1)- wf(end))./wf.^2);
t = linspace(-dt*(Np-1),dt*(Np-1),Np);

Ret = ifft(scan);
Ret = ifftshift(Ret);
phase = unwrap(angle(Ret));
Inten = abs(Ret).^2./max(abs(Ret).^2);



end

