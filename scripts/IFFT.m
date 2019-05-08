%% This function performs an inverse FFT of a e.g. time depending signal into the frequency space. The sampling rate is given
%  over the data interval 


function [xout,funcout]=IFFT(xvalues,yvalues)
M=length(xvalues);
Fsampling=1./((xvalues(M)-xvalues(1))./M);
funcout = ifft(ifftshift(yvalues)).*M;
xout=abs(Fsampling)./2*linspace(-1,1,M);
% xout=xout-xout(floor(NFFT/2)+1);

% % % netpo2 code
% % M=length(xvalues);
% % NFFT = 2.^nextpow2(M);
% % Fsampling=1./((xvalues(M)-xvalues(1))./M);
% % funcout = ifft(ifftshift(yvalues),NFFT).*M;
% % xout=abs(Fsampling)./2*linspace(-1,1,NFFT);

end
