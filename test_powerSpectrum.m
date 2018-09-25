%test my power spectrum vs toolbox version

dt=0.1;
t=0:dt:31.4;
y=sin(2*pi*w(1)*t)+sin(2*pi*w(2)*t);

% eix=1;
% ix=~(expt(eix).t>38);
% t=expt(eix).t(ix);
% y=expt(eix).R(ix,1);
% y=y-mean(y);

fs=1./dt;

[P,f]=pspectrum(y,fs,'Leakage',1);
% [P,f]=pwelch(y,[],[],4096,fs);
% [P,f]=periodogram(y,[],4096,fs);
[P1,f1]=powerSpectrum(y,fs);

figure(1); clf
plot(f,pow2db(P),f1,pow2db(P1));