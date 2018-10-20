function [acf,tau,peaklags]=autocorr(t,X,maxLag,peakThr,doPlot)
%autocorrelation function for each column of X, and finding peaks

%make vector column
if any(size(X)==1)
    X=X(:);
end

nX=size(X,2);
dt=mode(diff(t));

% if ~exist('params','var')||isempty(params)
%     params.maxLag=20;
%     params.peakThr=0.05;
% end
% maxLag=round(params.maxLag/dt);
% peakThr=params.peakThr;

if ~exist('maxLag','var')||isempty(maxLag)
    maxLag=t(end)/10;
end
if ~exist('peakThr','var')||isempty(peakThr)
    peakThr=0.05;
end
if ~exist('doPlot','var')
    doPlot=0;
end

maxLag=round(maxLag/dt);

%autocorrelation
acf=zeros(2*maxLag+1,nX);
lags=zeros(2*maxLag+1,nX);
for i=1:nX
    [acf(:,i),lags(:,i)]=xcorr(X(:,i),maxLag);
end

pts=peak_detector(lags*dt,acf,peakThr);

peaklags=zeros(1,nX);
for i=1:nX
    posix=pts(i).tMax>0;
%     posix=pts(i).tMax>0&pts(i).xMax>0;
    pospeaks=pts(i).tMax(posix);
    maxix=1; %first peak
%     [~,maxix]=max(pts(i).xMax(posix));
    peaklags(i)=pospeaks(maxix);
end

tau=lags*dt;
        
%plot to show result
if nargout==0 || doPlot==1

[P,f,pmax,fmax]=powerSpectrum(X,1./dt);

% [p,f]=pwelch(X,[],[],4096,1./dt);
% [~,ix]=max(p);
% fmax=f(ix);

Tpsd=1./fmax;

% ptsPSD=peak_detector(f,P,0);

tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        
        plot(tau(:,tix),acf(:,tix),'k','LineWidth',1);
        hold on
        plot(pts(tix).tMax,pts(tix).xMax,'rv')
%         plot(ptsPSD(tix).tMax,1./ptsPSD(tix).xMax,'bv')
        psdix=find(tau(:,tix)>=Tpsd(tix),1,'first');
        plot(Tpsd(tix),acf(psdix,tix),'bv')
        hold off  
        xlabel('lag')         
        ylabel('autocorrelation');
        axis tight
        
        YLIM=ylim();
        ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
        
    end

    function keypressFcn(~,event)
        switch(event.Key)
            case {'leftarrow'}
                if tix>1
                    tix=tix-1;
                    plotData()
                end
            case {'rightarrow'}
                if tix<nX
                    tix=tix+1;
                    plotData()
                end
        end
        
    end
    
end