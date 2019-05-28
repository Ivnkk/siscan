function intshg=SHGv2(dphase,int0,phase0)
siz=size(dphase);
amp0=sqrt(int0).*exp(1i*phase0);
n=length(phase0);
% l=1:2*n;
% id=mod(l,2)==1;
% m=2;
% fshg0=linspace(2*f0(1),2*f0(end),2*length(f0)-1);
% parfor i=1:siz(1)
% %     intt=ifft(amp0.*exp(1i*dphase(i,:)));
% %     shg0=intt.^2;
% %     shg=fft(shg0);
%     amp=amp0.*exp(1i*dphase(i,:));
%     shg=conv(amp,amp,'full');
%     int=abs(shg).^2;
% 	intshg(i,:)=sampling(fshg0,int,fshg,'');
% %     id=intshg(i,:)>=(max(intshg(i,:))/2+500);
% %     width(i)=sum(id)*(fshg(2)-fshg(1));
% end

% parfor i=1:siz(1)
%     intt=ifft(amp0.*exp(1i*dphase(i,:)),m*n);
%     shg0=intt.^2;
%     shg=fft(shg0);
%     int=abs(shg(1:n*m/2)).^2;
%     intshg(i,:)=sampling(fshg0,int,fshg,'');
% end

amp1=wextend('ar','sp0',amp0,siz(1)-1,'d');
intt=ifft(amp1.*exp(1i*dphase),2*n,2);
shg0=intt.^2;
shg=fft(shg0,[],2);
intshg=abs(shg(:,1:2:2*n)).^2;
intshg=intshg/max(max(intshg));