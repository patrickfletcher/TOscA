function [XNORM, XDT, XFILT, tfilt,XT]=preprocess(t,X,params,doPlot)
% convenience function with plotting to show whole pipeline

if ~exist('doPlot','var')
    doPlot=0;
end     

XNORM=normalizeTraces(t,X,params.norm.method,params.norm.methodpar);
        

[XDT,XT]=detrendTraces(t,XNORM,params.trend.method,params.trend.methodpar);

% if ~isfield(params.filt,'doTrim')
%     params.filt.doTrim=0;
% end
[XFILT,tfilt]=filterTraces(t,XDT,params.filt.method,params.filt.methodpar);

%plot to show result
if nargout==0 || doPlot==1

    fs=1/mode(diff(t));
    
    [PN,f]=pspectrum(XNORM,fs,'Leakage',0.9);
    [PDT,~]=pspectrum(XDT,fs,'Leakage',0.9);
    [PF,~]=pspectrum(XFILT,fs,'Leakage',0.9);
    
    PN=pow2db(PN);
    PDT=pow2db(PDT);
    PF=pow2db(PF);
    
    nX=size(X,2);
    tix=1;oldtix=1;
    figure('KeyPressFcn',@keypressFcn);
    
    ax(1)=subplot(4,1,1);
    hR=plot(t,X,'Color',0.75*[1,1,1]);
    grid on
    ylabel('raw')
    axis tight

    ax(2)=subplot(4,1,2);
    hN=plot(t,XNORM,'Color',0.75*[1,1,1]);
    hold on
    hT=plot(t,XT(:,tix),'r');
    hold off
    grid on
    ylabel('norm,trend')
    axis tight

    ax(3)=subplot(4,1,3);
    hDT=plot(t,XDT,'Color',0.75*[1,1,1]);
    hold on
    hF=plot(t,XFILT(:,tix),'r');
    hold off
    grid on
    xlabel('Time')
    ylabel('detrend,filt')
    axis tight

    ax(4)=subplot(4,1,4);
    hP=plot(f,PN(:,tix),f,PDT(:,tix),f,PF(:,tix));
    grid on
    xlabel('f')
    ylabel('power (dB)')
    axis tight
    legend('norm','detrend','filt')
    xlim([0,1])
    
    updateTrace()
    
end


%nested functions can see variables in caller's scope
    function updateTrace()
        
        %gray out old traces
        hR(oldtix).Color=0.75*[1,1,1];
        hN(oldtix).Color=0.75*[1,1,1];
        hDT(oldtix).Color=0.75*[1,1,1];
        hR(oldtix).ZData=zeros(size(t));
        hN(oldtix).ZData=zeros(size(t));
        hDT(oldtix).ZData=zeros(size(t));
        
        %black for new traces
        hR(tix).Color='k';
        hN(tix).Color='k';
        hDT(tix).Color='k';
        hR(tix).ZData=ones(size(t));
        hN(tix).ZData=ones(size(t));
        hDT(tix).ZData=ones(size(t));
        
        %trend and filter
        hT.YData=XT(:,tix);
        hF.YData=XFILT(:,tix);
        
        %power spectra
        hP(1).YData=PN(:,tix);
        hP(2).YData=PDT(:,tix);
        hP(3).YData=PF(:,tix);
        
        
    end

    function keypressFcn(~,event)
        switch(event.Key)
            case {'leftarrow'}
                if tix>1
                    oldtix=tix;
                    tix=tix-1;
                    updateTrace()
                end
            case {'rightarrow'}
                if tix<nX
                    oldtix=tix;
                    tix=tix+1;
                    updateTrace()
                end
        end
        
    end

end