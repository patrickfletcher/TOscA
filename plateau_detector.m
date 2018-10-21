function [points,features]=plateau_detector(t, X, varargin)
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
%  points - struct containing special points, for plotting

%check inputs and parse optional parameters

periodMarker='min';

if isvector(X)
    X=X(:);
end
    
[thrPtiles,fracUp,fracDown,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin{:});

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

%TODO: make sure initialization of upstate toggle is robust
DX=slopeY(t,X); %todo: noise-robust method; measure between tmin(i) and tmin(i+1)?
% isUp=ones(1,nX); %forces first point to be downward transition
% isUp=zeros(1,nX); %forces first point to be upward transition
isUp=X(1,:)>thrUp;
% isUp=X(1,:)>=thrUp | (X(1,:)<=thrUp&X(1,:)>=thrDown&DX(1,:)<0);

points(nX)=struct('tPer',[],'xPer',[],'iPer',[],...
    'tUp',[],'xUp',[],'iUp',[],'tDown',[],'xDown',[],'iDown',[],...
    'tMax',[],'xMax',[],'iMax',[],'tMin',[],'xMin',[],'iMin',[]);
%     'tDXMax',[],'xDXMax',[],'dxMax',[],'iDXMax',[],...
%     'tDXMin',[],'xDXMin',[],'dxMin',[],'iDXMin',[]);
for i=2:length(t)
    
    %if strated up, then went below thrUp but not below thrDown, set isUp=0
%     firstDown=isUp & X(i,:)<thrUp & arrayfun(@(x)isempty(x.iUp),points);
%     isUp(firstDown)=0;
    
    %check for downward-transition and upward-transition
    isDownTransition=isUp&X(i,:)<=thrDown;
    isUpTransition=~isUp&X(i,:)>=thrUp;
    
    if any(isDownTransition)
        idx=find(isDownTransition);
        for j=idx
            %index
            points(j).iDown(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrDown(j);
                thisT=(thrDown(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).tDown(end+1)=thisT;
            points(j).xDown(end+1)=thisX;
        end
        isUp(isDownTransition)=0;
    end
    
    if any(isUpTransition)
        idx=find(isUpTransition);
        for j=idx
            points(j).iUp(end+1)=i;
            
            if doInterp
                %interpolate t
                thisX=thrUp(j);
                thisT=(thrUp(j)-X(i-1,j))*(t(i)-t(i-1))/(X(i,j)-X(i-1,j))+t(i-1);
            else
                thisT=t(i);
                thisX=X(i,j);
            end
            
            points(j).tUp(end+1)=thisT;
            points(j).xUp(end+1)=thisX;
        end
        isUp(isUpTransition)=1;
    end
    
end

numUp=arrayfun(@(x) length(x.tUp),points);

% features=struct('period',[],'APD',[],'PF',[],'amp',[]);
features=struct();
for i=1:nX
        
    if numel(points(i).tUp)>1
    
%     trim
    if points(i).tDown(1)<=points(i).tUp(1)
        points(i).iDown=points(i).iDown(2:end);
        points(i).tDown=points(i).tDown(2:end);
        points(i).xDown=points(i).xDown(2:end);
    end
    if points(i).tDown(end)>=points(i).tUp(end)
        points(i).iDown=points(i).iDown(1:end-1);
        points(i).tDown=points(i).tDown(1:end-1);
        points(i).xDown=points(i).xDown(1:end-1);
    end
    
    points(i).iPer=points(i).iUp;
    points(i).tPer=points(i).tUp;
    points(i).xPer=points(i).xUp;
    
    nT=length(points(i).tUp)-1;
    
    features(i).range=globalXamp(i);
    features(i).thrUp=thrUp(i);
    features(i).thrDown=thrDown(i);
    features(i).period=diff(points(i).tUp);
    features(i).APD=points(i).tDown-points(i).tUp(1:end-1); %active phase duration
    features(i).PF=features(i).APD./features(i).period;
    
    for j=1:nT
        tt=points(i).iUp(j):points(i).iUp(j+1)-1;
        [xmax,imax]=max(X(tt,i));
        points(i).iMax(j)=tt(1)+imax-1;
        points(i).tMax(j)=t(tt(imax));
        points(i).xMax(j)=xmax;
        
        [xmin,imin]=min(X(tt,i));
        points(i).iMin(j)=tt(1)+imin-1;
        points(i).tMin(j)=t(tt(imin));
        points(i).xMin(j)=xmin;
        
%         [dxmax,idxmax]=max(DX(tt,i));
%         points(i).iDXMax(j)=tt(1)+idxmax-1;
%         points(i).tDXMax(j)=t(tt(idxmax));
%         points(i).xDXMax(j)=X(tt(idxmax),i);
%         points(i).dxMax(j)=dxmax;
%         
%         [dxmin,idxmin]=min(DX(tt,i));
%         points(i).iDXMin(j)=tt(1)+idxmin-1;
%         points(i).tDXMin(j)=t(tt(idxmin));
%         points(i).xDXMin(j)=X(tt(idxmin),i);
%         points(i).dxMin(j)=dxmin;
    end
    
    features(i).amp=points(i).xMax-points(i).xMin;
    
    else
        %had less than two minima: can't compute features.
        features(i).period=0;
        features(i).baseline=0;
        features(i).peaks=0;
        features(i).amp=0;
        features(i).APD=0;
        features(i).PF=0;
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
        axis tight
        hold on
        plot(points(tix).tUp,points(tix).xUp,'bs','Tag','plateau_detector')
        plot(points(tix).tDown,points(tix).xDown,'bo','Tag','plateau_detector')
        plot(points(tix).tMax,points(tix).xMax,'rv','Tag','plateau_detector')
        plot(points(tix).tMin,points(tix).xMin,'r^','Tag','plateau_detector')
%         plot(points(tix).tDXMax,points(tix).xDXMax,'g>','Tag','plateau_detector')
%         plot(points(tix).tDXMin,points(tix).xDXMin,'g<','Tag','plateau_detector')
        plot(xlim(),thrUp(tix)*[1,1],'r--','Tag','plateau_detector')
        if thrUp~=thrDown
            plot(xlim(),thrDown(tix)*[1,1],'b--','Tag','plateau_detector')
        end
        xlabel('t')
        ylabel('x')
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

function [thrPtiles,fracUp,fracDown,doInterp,doPlot,figID,dokeypress]=parseArgs(t, X, varargin)

    %default parameters
    defaultF=[0.5,0.4]; %near halfmax
    thrPtiles=[0,100];
    doInterp=true;
    doPlot=false;
    doKeypress=true;
    figID=[];
    
    p=inputParser;
    validX=@(x) isreal(x) && size(x,1)==length(t); %traces are columns of X
    validSwitch=@(x) isscalar(x) && (isnumeric(x)||islogical(x));
    addRequired(p,'t',@(x) isreal(x));
    addRequired(p,'X',validX);
    addOptional(p,'f',defaultF);  %somehow adding the validation function here messes things up
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