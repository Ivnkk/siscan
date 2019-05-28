function shg=SHGv22(dphase,amp0)
siz=size(dphase);
amp1=wextend('ar','sp0',amp0,siz(1)-1,'d');
intt=ifft(amp1.*exp(1i*dphase),[],2);
shg0=intt.^2;
shg=fft(shg0,[],2);
