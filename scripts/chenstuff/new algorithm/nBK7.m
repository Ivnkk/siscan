function n=nBK7(wl)
%Index of refraction of BK7 as function of wavelength in \mum
n2=1+1.03961212*wl.^2./(wl.^2-0.00600069867)+0.231792344*wl.^2./(wl.^2-0.0200179144)+1.01046945*wl.^2./(wl.^2-103.560653);
id=n2>=0.1;
n2=id.*n2;
n=sqrt(n2);