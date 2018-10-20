function Xbar=weightedTraceMean(t,X,doPlot)
% "shrinkage" mean idea: weight each column of X by its inverse distance to the mean trace.

nX=size(X,2);

if ~exist('doPlot','var')
    doPlot=0;
end

Xbar1=mean(X,2);
XS=zeros(size(X));
shifts=zeros(1,nX);
for i=1:nX
    D=norm(X(:,i)-Xbar1, 1);
    w(i)=1/D;
end
Xbar=sum(X.*w,2)/sum(w);

        

%plot to show result
if nargout==0 || doPlot==1

tix=1;
figure('KeyPressFcn',@keypressFcn);
plotData()
    
end


%nested functions can see variables in caller's scope
    function plotData()
        
%         subplot(2,1,1)
        plot(t,Xbar1,'r')
        hold on
        plot(t,Xbar,'k')
        hold off
        
        xlabel('t')         
        ylabel('X');
        axis tight
        YLIM=ylim();
        ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
        
        
%         subplot(2,1,2)
%         plot(t,XS,'color',[0.6,0.6,0.6])
%         hold on
%         plot(t,XS(:,tix),'b')
%         plot(t,Xbar,'k','linewidth',2)
%         hold off
%         
%         xlabel('t')         
%         ylabel('X');
%         axis tight
%         YLIM=ylim();
%         ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
        
        
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