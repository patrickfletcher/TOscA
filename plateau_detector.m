function [trace,features]=plateau_detector(t, X, varargin)
%PLATEAU_DETECTOR A plateau detector for oscillating timeseries. Uses Shmitt
% trigger concept - upward transition threshold >= downward transition
% threshold. The timeseries is divided into periods measured from upward
% transition to upward transition. 
%
% If results are plotted, use left/right arrows to switch between columns
% of X.
%
% Required Inputs:
%  t - vector of timepoints
%  X - vector matching t, or matrix with columns matching t
%
% Optional input:
%  f - fraction of amplitude for threshold: 
%     scalar {default, f=0.5}, upward=downward threshold,
%     vector, f=[upFraction,downFraction] with upFraction>=downFraction.
%
% optional name-value pair arguments:
%  'Interpolate', {true}/false
%  'ThresholdPercentiles', default=[0,100] (threshold is fraction of this range of amplitudes)
%  'Plot',true/{false} (use to control plotting with output arguments)
%
% Outputs:
%  features - struct containg feature distributions for each column of X
%  trace - struct containing special points, for plotting

%check inputs and parse optional parameters
[thrPtiles,fracUp,fracDown,doInterp,doPlot]=parseArgs(t, X, varargin{:});

nX=size(X,2); %number of traces

globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

%independent thresholds for each trace
thrUp=globalXmin + fracUp*globalXamp;
if fracDown>0
    thrDown=globalXmin + fracDown*globalXamp;
else
    thrDown=thrUp;
end

isUp=X(1,:)>thrUp;

DX=slopeY(t,X); %todo: noise-robust method; measure between tmin(i) and tmin(i+1)?

trace(nX)=struct('tUp',[],'xUp',[],'iUp',[],'tDown',[],'xDown',[],'iDown',[],...
                 'tMax',[],'xMax',[],'iMax',[],'tMin',[],'xMin',[],'iMin',[],...
                 'tDXMax',[],'xDXMax',[],'dxMax',[],'iDXMax',[],...
                 'tDXMin',[],'xDXMin',[],'dxMin',[],'iDXMin',[]);
for i=2:length(t)
    
    %check for downward-transition and upward-transition
    isDownTransition=isUp&X(i,:)<thrDown;
    isUpTransition=~isUp&X(i,:)>thrUp;
    
    if any(isDownTransition)
        idx=find(isDownTransition);
        for j=idx
            %index
            trace(j).iDown=[trace(j).iDown,i];
            
            if doInterp
                %interpolate t
                thisX=thrDown(j);
                thisT=(thrDown(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            trace(j).tDown=[trace(j).tDown,thisT];
            trace(j).xDown=[trace(j).xDown,thisX];
        end
        isUp(isDownTransition)=0;
    end
    
    if any(isUpTransition)
        idx=find(isUpTransition);
        for j=idx
            trace(j).iUp=[trace(j).iUp,i];
            
            if doInterp
                %interpolate t
                thisX=thrUp(j);
                thisT=(thrUp(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            trace(j).tUp=[trace(j).tUp,thisT];
            trace(j).xUp=[trace(j).xUp,thisX];
        end
        isUp(isUpTransition)=1;
    end
    
end

%remove downward transitions that are outside the first or last upward
%transition
firstIsDown=arrayfun(@(x) x.tDown(1)<x.tUp(1),trace);
lastIsDown=arrayfun(@(x) x.tDown(end)>x.tUp(end),trace);

idx=find(firstIsDown);
for j=idx
    trace(j).iDown=trace(j).iDown(2:end);
    trace(j).tDown=trace(j).tDown(2:end);
    trace(j).xDown=trace(j).xDown(2:end);
end
idx=find(lastIsDown);
for j=idx
    trace(j).iDown=trace(j).iDown(1:end-1);
    trace(j).tDown=trace(j).tDown(1:end-1);
    trace(j).xDown=trace(j).xDown(1:end-1);
end

%append trace dependent info?

%features
for i=1:nX
    nT=length(trace(i).tUp)-1;
    features(i).T=diff(trace(i).tUp);
    features(i).APD=trace(i).tDown-trace(i).tUp(1:end-1); %active phase duration
    features(i).PF=features(i).APD/features(i).T;
    for j=1:nT
        tt=trace(i).iUp(j):trace(i).iUp(j+1)-1;
        [xmax,imax]=max(X(tt,i));
        trace(i).iMax(j)=imax;
        trace(i).tMax(j)=t(tt(imax));
        trace(i).xMax(j)=xmax;
        
        [xmin,imin]=min(X(tt,i));
        trace(i).iMin(j)=imin;
        trace(i).tMin(j)=t(tt(imin));
        trace(i).xMin(j)=xmin;
        
        [dxmax,idxmax]=max(DX(tt,i));
        trace(i).iDXMax(j)=idxmax;
        trace(i).tDXMax(j)=t(tt(idxmax));
        trace(i).xDXMax(j)=X(tt(idxmax),i);
        trace(i).dxMax(j)=dxmax;
        
        [dxmin,idxmin]=min(DX(tt,i));
        trace(i).iDXMin(j)=idxmin;
        trace(i).tDXMin(j)=t(tt(idxmin));
        trace(i).xDXMin(j)=X(tt(idxmin),i);
        trace(i).dxMin(j)=dxmin;
        
        features(i).amp(j)=xmax-xmin;
    end
end


%plot to show performance
if nargout==0||doPlot==1
    tix=1;
    figure('KeyPressFcn',@keypressFcn);
    plotData()
end

%nested functions can see variables in caller's scope
    function plotData()
        clf
        plot(t,X(:,tix),'k-')
        hold on
        plot(trace(tix).tUp,trace(tix).xUp,'bs')
        plot(trace(tix).tDown,trace(tix).xDown,'bo')
        plot(trace(tix).tMax,trace(tix).xMax,'rv')
        plot(trace(tix).tMin,trace(tix).xMin,'r^')
%         plot(trace(tix).tDXMax,trace(tix).xDXMax,'g>')
%         plot(trace(tix).tDXMin,trace(tix).xDXMin,'g<')
        plot(xlim(),thrUp(tix)*[1,1],'r--')
        if thrUp~=thrDown
            plot(xlim(),thrDown(tix)*[1,1],'b--')
        end
        xlabel('t')
        ylabel('x')
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

function [thrPtiles,fracUp,fracDown,doInterp,doPlot]=parseArgs(t, X, varargin)

    %default parameters
    defaultF=0.5; %near halfmax
    thrPtiles=[0,100];
    doInterp=true;
    doPlot=false;
    
    p=inputParser;
    validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
    validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
    addRequired(p,'t',@(x) isreal(x));
    addRequired(p,'X',validX);
    addOptional(p,'f',defaultF);  %somehow adding the validation function here messes things up
    addParameter(p,'ThresholdPercentiles',thrPtiles,@(x) isreal(x) && numel(x)==2);
    addParameter(p,'Interpolate',doInterp,validSwitch);
    addParameter(p,'Plot',doPlot,validSwitch);
    
    parse(p,t,X,varargin{:});
    
    f=p.Results.f;
    thrPtiles=p.Results.ThresholdPercentiles;
    doInterp=p.Results.Interpolate;
    doPlot=p.Results.Plot;
    
    if isscalar(f)
        fracUp=f;
        fracDown=f;
    elseif length(f)==2
        fracUp=f(1);
        fracDown=f(2);
    else
        error('Invalid value for fractional threshold')
    end
    
    %snap to fracUp?
    if fracDown>fracUp
        fracDown=fracUp;
    end
    
end