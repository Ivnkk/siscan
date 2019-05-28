function intshg=SHGv21(dphase,int0,phase0)
siz=size(dphase);
amp0=sqrt(int0).*exp(1i*phase0);
amp1=wextend('ar','sp0',amp0,siz(1)-1,'d');
intt=ifft(amp1.*exp(1i*dphase),[],2);
shg0=intt.^2;
shg=fft(shg0,[],2);
intshg=abs(shg).^2;