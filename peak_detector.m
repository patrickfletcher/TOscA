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

delta=0.5;
deltaIsFraction=true;
thrPtiles=[0,100];
activeFraction=0.5;

nX=size(X,2); %number of traces

globalXmin=prctile(X,thrPtiles(1),1);
globalXmax=prctile(X,thrPtiles(2),1);
globalXamp=globalXmax-globalXmin;

if deltaIsFraction
    delta=delta*globalXamp;
end

hold off; plot(t,X); hold on

points(nX)=struct('range',[],...
    'tMax',[],'xMax',[],'iMax',[],'tMin',[],'xMin',[],'iMin',[]);

maxIsNext=zeros(1,nX); %always look for a minimum first
lastMax=X(1,:); %initialize first point
lastMin=X(1,:);
maxix=ones(1,nX);
minix=ones(1,nX);
for i=2:length(t)
    this=X(i,:);
    [lastMax,maxrow]=max([this;lastMax],[],1); %update most recent maximum for all traces
    maxix(maxrow==1)=i; %maxrow==1 means new max was found, maxrow==2 if lastMax was higher.
    [lastMin,minrow]=min([this;lastMin],[],1);
    minix(minrow==1)=i;
    
    maxFound=maxIsNext & this<lastMax-delta;
    idx=find(maxFound); %for each trace with max found, store the point
    for j=idx
        points(j).iMax(end+1)=maxix(j);
        points(j).tMax(end+1)=t(maxix(j));
        points(j).xMax(end+1)=lastMax(j);
        lastMin(maxFound)=this;
        minix(maxFound)=i;
%         lastMax(maxFound)=-Inf;
        maxIsNext(maxFound)=false;
        plot(points(j).tMax(end),points(j).xMax(end),'v');
    end
    
    minFound=~maxIsNext & this>lastMin+delta; 
    idx=find(minFound);
    for j=idx
        points(j).iMin(end+1)=minix(j);
        points(j).tMin(end+1)=t(minix(j));
        points(j).xMin(end+1)=lastMin(j);
%         lastMin(minFound)=Inf;
        lastMax(minFound)=this;
        maxix(minFound)=i;
        maxIsNext(minFound)=true;
        plot(points(j).tMin(end),points(j).xMin(end),'^');
    end
    
end

for j=1:nX
plot(points(j).tMax,points(j).xMax,'v');
plot(points(j).tMin,points(j).xMin,'^');
end


features=struct('T',[],'APD',[],'PF',[],'amp',[]);
for i=1:nX
    points(i).range=globalXamp(i);
    
    if numel(points(i).tMin)>1
        
    %trim leading and trailing maxima
    if points(i).tMax(1)<points(i).tMin(1)
        points(i).iMax=points(i).iMax(2:end);
        points(i).tMax=points(i).tMax(2:end);
        points(i).xMax=points(i).xMax(2:end);
    end
    if points(i).tMax(end)>points(i).tMin(end)
        points(i).iMax=points(i).iMax(1:end-1);
        points(i).tMax=points(i).tMax(1:end-1);
        points(i).xMax=points(i).xMax(1:end-1);
    end
    
    features(i).T=diff(points(i).tMin); %period defined from minimum to minimum
    
    %interpolate active phase threshold crossing per period
    %local linear detrend using two minima?
    nPer=length(trace(i).tUp)-1;
    for j=1:nPer
        
%         features(i).APD(j)=
    end
    
    features(i).PF=features(i).APD./features(i).T;
    features(i).amp=points(i).xMax-points(i).xMin;
   
    end
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