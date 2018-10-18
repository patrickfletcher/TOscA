function [points,features]=peak_detector(t, X, varargin)
%PEAK_DETECTOR A peak detector for oscillating timeseries. The timeseries is divided into periods measured between
%successive minima. Amplitude and time window sizes can be set to select scale of peaks.
%
% If results are plotted, use left/right arrows to switch between columns
% of X.
%
% Required Inputs:
%  t - vector of timepoints
%  X - vector matching t, or matrix with columns matching t
%
% Optional input:
%  delta - fraction of amplitude for threshold: 
%     - x is max if x>last min + delta
%     scalar {default, f=0.5}
%  wsz - window size within which to check for extrema
%
% Optional name-value pair arguments:
%  'Interpolate', {true}/false - for sub-resolution peak locations
%  'ThresholdPercentiles', default=[0,100] (threshold is fraction of this range of amplitudes)
%  'Plot',true/{false} (use to control plotting with output arguments)
%
% Outputs:
%  features - struct containg feature distributions for each column of X
%  trace - struct containing special points, for plotting

%check inputs and parse optional parameters
% [thrPtiles,delta,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

delta=0.25;
deltaIsFraction=true;
thrPtiles=[0,100];

nX=size(X,2); %number of traces

globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

if deltaIsFraction
    delta=delta*globalXamp;
end


points(nX)=struct('range',[],...
    'tMax',[],'xMax',[],'iMax',[],'tMin',[],'xMin',[],'iMin',[]);

maxIsNext=zeros(1,nX); %always look for a minimum first
lastMax=X(1,:); %initialize first point
lastMin=X(1,:);
for i=2:length(t)
    this=X(i,:);
    lastMax=max(this,lastMax);
    lastMin=min(this,lastMin);
    
    maxCondition=maxIsNext & this<lastMax-delta;
    idx=find(maxCondition);
    for j=idx
        points(j).iMax(end+1)=i;
        points(j).tMax(end+1)=t(i);
        points(j).xMax(end+1)=lastMax(j);
    end
    lastMin(maxCondition)=this(maxCondition);
    
    minCondition=~maxIsNext & this>lastMin+delta; 
    idx=find(minCondition);
    for j=idx
        points(j).iMin(end+1)=i;
        points(j).tMin(end+1)=t(i);
        points(j).xMin(end+1)=lastMin(j);
    end
    lastMax(minCondition)=this(minCondition);
    
end



features=struct('T',[],'APD',[],'PF',[],'amp',[]);
for i=1:nX
    
    points(i).range=globalXamp(i);
    points(i).thrUp=thrUp(i);
    points(i).thrDown=thrDown(i);
    
    if numel(points(i).tUp)>1
        
    %trim so tups contain all tdowns
    if points(i).tDown(1)<points(i).tUp(1)
        points(i).iDown=points(i).iDown(2:end);
        points(i).tDown=points(i).tDown(2:end);
        points(i).xDown=points(i).xDown(2:end);
    end
    if points(i).tDown(end)>points(i).tUp(end)
        points(i).iDown=points(i).iDown(1:end-1);
        points(i).tDown=points(i).tDown(1:end-1);
        points(i).xDown=points(i).xDown(1:end-1);
    end
    
    nT=length(points(i).tUp)-1;
    features(i).T=diff(points(i).tUp);
    features(i).APD=points(i).tDown-points(i).tUp(1:end-1); %active phase duration
    features(i).PF=features(i).APD./features(i).T;
   
    end
    features(i).amp=points(i).xMax-points(i).xMin;
end


%plot to show performance
if nargout==0||doPlot==1
    tix=1;
    if isempty(figID)
        figID=gcf; %new figure if none available, otherwise current fig
    else
        figID=figure(figID);
    end
    if dokeypress
        figID.KeyPressFcn=@keypressFcn;
    end
    
    plotData()
end

%nested functions can see variables in caller's scope
    function plotData()
%         clf
        delete(findobj(gca,'Tag','plateau_detector'));

        plot(t,X(:,tix),'k-','Tag','plateau_detector')
        hold on
        plot(points(tix).tUp,points(tix).xUp,'bs','Tag','plateau_detector')
        plot(points(tix).tDown,points(tix).xDown,'bo','Tag','plateau_detector')
        plot(points(tix).tMax,points(tix).xMax,'rv','Tag','plateau_detector')
        plot(points(tix).tMin,points(tix).xMin,'r^','Tag','plateau_detector')
%         plot(trace(tix).tDXMax,trace(tix).xDXMax,'g>','Tag','plateau_detector')
%         plot(trace(tix).tDXMin,trace(tix).xDXMin,'g<','Tag','plateau_detector')
        plot(xlim(),thrUp(tix)*[1,1],'r--','Tag','plateau_detector')
        if thrUp~=thrDown
            plot(xlim(),thrDown(tix)*[1,1],'b--','Tag','plateau_detector')
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

function [thrPtiles,fracUp,fracDown,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

    %default parameters
    defaultDelta=0.25;
    doInterp=true;
    doPlot=false;
    doKeypress=false;
    figID=[];
    
    p=inputParser;
    validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
    validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
    addRequired(p,'t',@(x) isreal(x));
    addRequired(p,'X',validX);
    addOptional(p,'delta',defaultDelta);
    addParameter(p,'ThresholdPercentiles',thrPtiles,@(x) isreal(x) && numel(x)==2);
    addParameter(p,'Interpolate',doInterp,validSwitch);
    addParameter(p,'Plot',doPlot,validSwitch);
    addParameter(p,'Keypress',doKeypress,validSwitch);
    addParameter(p,'FigureID',figID);
    
    parse(p,t,X,varargin{:});
    
    f=p.Results.f;
    thrPtiles=p.Results.ThresholdPercentiles;
    doInterp=p.Results.Interpolate;
    doPlot=p.Results.Plot;
    figID=p.Results.FigureID;
    dokeypress=p.Results.Keypress;
    
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