function [Xbar,shifts]=alignToMean(t,X,maxLag,doPlot)
%iterative timeshifting of columns of X toward the mean
% "shrinkage" averaging of timeseries to remove small phase shifts?

nX=size(X,2);

if ~exist('maxLag','var')||isempty(maxLag)
    maxLag=t(end)/10;
end
if ~exist('doPlot','var')
    doPlot=0;
end

dt=mode(diff(t));
maxLag=round(maxLag/dt);

Xbar1=mean(X,2);
XS=zeros(size(X));
shifts=zeros(1,nX);
for i=1:nX
    [xa,ya,D]=alignsignals(Xbar1,X(:,i),maxLag);
    shifts(i)=D;
    if D>0
        XS(:,i)=[ya(D+1:end);zeros(D,1)];
    elseif D==0
        XS(:,i)=X(:,i); %no change
    else
        D=-D;
        XS(:,i)=ya(1:end-D);
    end
end
% Xbar=mean(XS,2);
Xbar=weightedTraceMean(t,XS);
        

%plot to show result
if nargout==0 || doPlot==1

tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        
        subplot(2,1,1)
        plot(t,X,'color',[0.6,0.6,0.6])
        hold on
        plot(t,X(:,tix),'b')
        plot(t,Xbar1,'k','linewidth',2)
        hold off
        
        xlabel('t')         
        ylabel('X');
        axis tight
        YLIM=ylim();
        ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
        
        
        subplot(2,1,2)
        plot(t,XS,'color',[0.6,0.6,0.6])
        hold on
        plot(t,XS(:,tix),'b')
        plot(t,Xbar,'k','linewidth',2)
        hold off
        
        xlabel('t')         
        ylabel('X');
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