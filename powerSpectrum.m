function [P1,f,Pmax,fmax]=powerSpectrum(x,fs,params)
%if x=matrix, columns are timeseries

%actually this is a periodogram

%make vector column
if any(size(x)==1)
    x=x(:);
end

doplot=0;
if exist('figID','var')&&~isempty(figID)
    doplot=1;
end

if ~exist('params','var')||isempty(params)
    params.n=4096;
end
n=params.n;

L=size(x,1); %assume column is timeseries

f=fs*(0:n/2-1)/n;

Y=fft(x,n);

% phaseY = unwrap(angle(Y));

%Power of each freq components 
% P2=Y.*conj(Y)/(n*L);
P2=Y.*conj(Y)/(fs*L);
% P2=abs(Y).^2/(n*L);
P1=P2(1:n/2,:);
P1(2:end-1,:)=2*P1(2:end-1,:);

[Pmax,ix]=max(P1);

fmax=f(ix);

if nargout==0
figure(figID);
plot(f,pow2db(P1),'LineWidth',1); 
% set(gca,'yscale','log')
xlim([0,1.5])
title('One Sided Power Spectral Density');       
xlabel('frequency')         
ylabel('power');

end