%% This function performs a FFT of a e.g. time depending signal into the frequency space. The sampling rate is given
%  over the data interval 

function [xout,funcout]=FFT(xvalues,yvalues)
M=length(xvalues);
Fsampling=1/((xvalues(M)-xvalues(1))./M);
funcout = fft(yvalues);
xout = abs(Fsampling).*linspace(0,1,M);
end
