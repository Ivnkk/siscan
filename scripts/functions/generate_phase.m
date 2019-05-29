function [f_axis,phase] = generate_phase(wl,wl_center,phase_vec)
%GENERATE_PHASE Generate phase vector with given phase
% wl and wl_center in nm, phase in units of fs^n as polyval vector.
% Example: generate_phase(600:900,800,[20,0,pi])
%          here a gdd of 20fs^2 is assumed and a CEP of pi
% Jan Vogelsang, August 2018


f_axis = 300./wl; %freq, rad/fs
f_center = 300./wl_center;
f_axis = linspace(f_axis(1),f_axis(end),length(f_axis)); %reinterpolate to have equal spacing on frequency axis

taylor_factor = 1./factorial(length(phase_vec)-1:-1:0);

phase = 2*pi*polyval(taylor_factor.*phase_vec,f_axis-f_center);

end

