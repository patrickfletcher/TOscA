classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    %TODO: plotting efficiently - plot all timeseries onces, and use an
    %update function to set the visibility on/off (also would support
    %option to plot out-of-focus traces as light gray)
    
    properties
        name=''
        date=''
        sex=''
        condition='' %eg. WT vs KO
        filename=''
        fullfile=''
        notes={}
        
        t
        dt
        
        %segment is a subinterval of the trace, the basic unit within which
        %periodicity will be detected and features will be computed
        segment=struct('name','','endpoints',[],'ix',[],...
            'points',struct(),'features_periods',struct(),'features_trace',struct())
        nS
        
        %TODO: support 3D array - 3rd dim is observable id (one [nT x nX] page per observable)
        % this requires supporting different pipeline options for each
        % observable and allows computation of new higher order features,
        % such as phase relationship
        % --> actually, will likely determine period markers using only one
        % observable, then compute features relative to those markers in
        % other traces?
        
        %TODO: Customizeable pipeline? eg. array processingStep objects to define
        %steps (name+operation+method+params), and pointer to which one
        %gets periods measured and which one features get measured. pages
        %of X correspond to each step: eg. X(:,:,1)=raw, X(:,:,2)=norm, ..., X(:,:,end)=per
        X
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        
        include %set element to zero to exclude a trace; X(:,include)
        nX
        
        featureFcn
        fnames_periods
        fnames_trace
        
        fs %1/dt
        f
        psd
        
        group
        nG
        Xg %average traces with same group (always use Xdetrend? probably)
        Xgfilt %filter using filterMethod, filterParam
        
        resultsTrace=table
        resultsPeriod=table
        
    end
    
    %store parameters used
    properties
        
%         tOmit
%         nOmit

%support variable length parameter lists, eg. for name-value pairs
        normMethod='none'
        normParam=[]
        trendMethod='none'
        trendParam=[]
        filterMethod='none'
        filterParam=[]
        
        averageMethod='arithmetic'
        averageParam=[]
        
        featureMethod='threshold'
        featureParam={}
        
        interpMethod='linear'
        
%         thrFrac=[0.55,0.45]
%         delta=0.25
%         Tbig=15
%         Tsmall=4
        
        %keep track of which trace is in focus for plotting - no, this should be a property of the figure
%         tix=1;

        fig_handles=matlab.ui.Figure.empty %array of figure handles spawned by this object - use to synchronize trace in focus?
        tix %in focus trace
        active_fig %index into fig_handles to maintain its focus upon updates
    end
    
    methods
        %constructor
        % Experiment(): requests file
        % Experiment(filename): looks for given file
        % Experiment(data): uses numeric matrix data=[t,X]
        % TODO: Experiment(___,Name,Value) : name-value pairs to set parameters, metadata, etc.
        function expt=Experiment(varargin)
            
            %parse inputs
            loadFile=true;
            fullfile=[];
            filename=[];
            if ~isempty(varargin) && ( ischar(varargin{1}) || isstring(varargin{1}) )
                [path,filename,ext]=fileparts(varargin{1});
                fullfile=[path,filesep,filename,ext];
            elseif ~isempty(varargin) && isnumeric(varargin{1})
                data=varargin{1};
                loadFile=false;
            end
            
            %load a new file
            if loadFile
                if isempty(fullfile) || ~exist(fullfile,'file')
                    [filename,path]=uigetfile({'*.xls*'});
                    fullfile=[path,filename];
                end
                
                [num,txt,raw]=xlsread(fullfile);
                %extract data
                ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
                data=cell2mat(raw(ixs:end,1:size(num,2)));
                
                %look for special metadata entries
                headerNames=txt(:,1);
                keyword={'name','date','sex','condition'};
                for i=1:length(keyword)
                    r=find(strcmpi(headerNames,keyword{i}));
                    if ~isempty(r)
                        expt.(keyword{i})=txt{r,2};
                    end
                end
                
                r=find(strcmpi(headerNames,'group'));
                if ~isempty(r)
                    expt.setGroup(cell2mat(raw(r,2:size(num,2))))
                end
                
            end
            
            expt.filename=filename;
            expt.fullfile=fullfile;
                
            time=data(:,1);
            time=time-time(1); %shift time to start at zero
            Xraw=data(:,2:end);
            
            DT=mode(diff(time)); 

            %interpolate any missing data using neighboring time points
            %TODO: make this optional?
            
            timeIsNan=isnan(time); %this would be a problem!
            if any(timeIsNan)
                expt.notes=[expt.notes,{'missing value in time column'}];
            end
            
            skip=any(abs(diff(time)-DT)>2*DT);
            if any(skip)
                expt.notes=[expt.notes,{'large timestep detected'}];
            end
            
            rowIsNan=any(isnan(Xraw),2);
            if any(rowIsNan)
                expt.notes=[expt.notes,{'some rows have NaN'}];
                r=find(rowIsNan);

                if r(1)==1
                    time=time(2:end,1);
                    Xraw=Xraw(2:end,:);
                    r=r(2:end)-1;
                end
                for ii=1:length(r)
                    for jj=1:size(Xraw,2)
                        iix=[r(ii)-1,r(ii)+1];
                        x=time(iix,1);
                        v=Xraw(iix,jj);
                        xq=time(r(ii),1);
                        Xraw(r(ii),jj)=interp1(x,v,xq,expt.interpMethod);
                    end
                end
            end
%             expt.t=time;
%             expt.X=Xraw;
            
            %interpolate on a fixed grid with given DT
            %TODO: DT as parameter (for same DT across expts)
            % - alt, make this a method with DT parameter
            XI=[];
            ti=0:DT:time(end); ti=ti';
            for j=1:size(Xraw,2)
                XI(:,j)=interp1(time,Xraw(:,j),ti,'pchip');
            end
            expt.t=ti;
            expt.X=XI;
            
            expt.nX=size(expt.X,2);
            
            expt.dt=DT;
            expt.fs=1/expt.dt;
            
            %default segment = whole trace
            expt.defineSegments({'1'},[ti(1),ti(end)]);
%             expt.setGroup(ones(1,expt.nX));
            expt.include=true(1,expt.nX);
            
            expt.tix=1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Setup helper functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function defineSegments(expt,names,startTimes)
            %expects names as cell array, times as vector
            expt.nS=length(names);
%             expt.segment(expt.nS).name=''; %expand the array of struct
            startTimes(end+1)=expt.t(end); %append final time
            for i=1:expt.nS
                expt.segment(i).name=names{i};
                expt.segment(i).endpoints=[startTimes(i),startTimes(i+1)];
                expt.segment(i).ix=expt.t>=startTimes(i) & expt.t<=startTimes(i+1);
            end
        end
        
        function setGroup(expt,groupvar)
            expt.group=groupvar;
            expt.nG=unique(groupvar);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO: forced order, but flags to choose which steps (if any) to apply?
        %     - alt: arbitrary list of operations, with fully qualified options.
        %TODO: per-interval, or whole trace options (eg. linear detrend per interval)
        %TODO: single preprocessing convenience function. how to handle all the options+params? 
        
        function normalize(expt,method,methodPar,doPlot)
            
            if ~exist('doPlot','var')
                doPlot=false;
            end
            expt.Xnorm=normalizeTraces(expt.t,expt.X,method,methodPar,doPlot);
            expt.normMethod=method;
            expt.normParam=methodPar;
        end
        
        %detrend supports per-interval: for piecewise linear detrend.
        function detrend(expt,method,methodPar,perSegment,flattenMethod,doPlot)
            if ~exist('flattenMethod','var')
                flattenMethod='none';
            end
            if ~exist('doPlot','var')
                doPlot=false;
            end
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
                expt.normMethod='none';
                expt.normParam=[];
            end
            if perSegment
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    [Xdt,Xt]=detrendTraces(expt.t(thisIx),expt.Xnorm(thisIx,:),method,methodPar,doPlot);
                    
                    if ~(flattenMethod=="none")
                        %apply flattening. TODO: same for all segments? stitch
                        %last values??
                        if i==1 
                            switch flattenMethod
                                case 'mean'
                                    addVal=mean(expt.Xnorm(thisIx,:),1);
                                case 'median'
                                    addVal=median(expt.Xnorm(thisIx,:),1);
                                case 'first'
                                    addVal=expt.Xnorm(find(thisIx,1,'first'),:);
                                case 'meanTrend'
                                    addVal=mean(Xt,1);
                                case 'medianTrend'
                                    addVal=median(Xt,1);
                                case 'firstTrend'
                                    addVal=Xt(1,:);
                                otherwise
                                    error(['Unknown flatten method: ', flattenMethod])
                            end
                        else
                            addVal=last-Xdt(1,:); %makes the trace continuous
                        end
                        Xdt=Xdt+addVal;
                        last=Xdt(end,:);
                    end
                    
                    expt.Xtrend(thisIx,:)=Xt;
                    expt.Xdetrend(thisIx,:)=Xdt;
                end
                
                %TODO: this plots a figure for each segment. could suppress & make combined plot? Really need a full
                %preprocessing helper function/plot...?
%                 expt.plotTrace('detrend');
                
            else
                [expt.Xdetrend,expt.Xtrend]=detrendTraces(expt.t,expt.Xnorm,method,methodPar,doPlot);
            end
            expt.trendMethod=method;
            expt.trendParam=methodPar;
        end
        
        function averageTraces(expt)
            XDTg=zeros(length(expt.t),expt.nG);
            for j=1:nI  
                XDTg(:,j)=mean(XDTg(:,I==uI(j)),2);  
            end
            expt.Xg=XDTg;
        end
        
        
        function filter(expt,method,methodPar,doPlot)
            if ~exist('doPlot','var')
                doPlot=false;
            end
            if isempty(expt.Xnorm)
                expt.Xnorm=expt.X;
                expt.normMethod='none';
                expt.normParam=[];
            end
            if isempty(expt.Xdetrend)
                expt.Xdetrend=expt.X;
                expt.trendMethod='none';
                expt.trendParam=[];
            end
            expt.Xfilt=filterTraces(expt.t,expt.Xdetrend,method,methodPar,doPlot);
            expt.filterMethod=method;
            expt.filterParam=methodPar;
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO expand this - detect_periods return only points/fcns, then compute_features should apply fcn to relevant
        %traces.
        function compute_features(expt,method,varargin)
            expt.detect_periods(method,varargin{:});
            expt.periodogram();
            expt.fnames_periods=fieldnames(expt.segment(1).features_periods)';
            expt.fnames_trace=fieldnames(expt.segment(1).features_trace)';
            expt.featureMethod=method;
            expt.featureParam=varargin; %methodpar should be varargin?
            
            expt.buildResultsTable();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plotTrace(expt,whichPlot,tix,showPts,doInteractive)
            %dispatch function for plotting expt data
            
            %TODO: support cell array for whichplot
                   
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('showPts','var')||isempty(showPts)
                showPts=true;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
            
            if ~iscell(whichPlot)
                whichPlot={whichPlot};
            end
            
            %interactive figure: set up callbacks
            if doInteractive
                figID=gcf;
                figID.Name=['Traces: ',expt.filename];
                figID.NumberTitle='off';
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
                figID.KeyPressFcn=@expt.commonKeypress;
                figID.UserData={'trace',whichPlot,showPts};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            nPlots=length(whichPlot);
            for i=1:nPlots
                subplot(nPlots,1,i)
                expt.plot_t(whichPlot{i},showPts)
            end
        end
        
        function plotPeriodogram(expt,tix,doInteractive)
            %overlay psd for each segment
            if ~isempty(expt.psd)
            
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
                      
            if doInteractive
                figID=gcf;
                figID.Name=['Periodogram: ',expt.filename];
                figID.NumberTitle='off';
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
                figID.KeyPressFcn=@expt.commonKeypress;
                figID.UserData={'psd'};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            expt.plot_psd();
            
            end
        end
        
        function plotFeatures(expt,xname,yname,tix,doInteractive)
            %TODO: how to manage plotting options?
                
%             if ~isempty(expt.segment(1).features)
            
            if exist('tix','var')&&~isempty(tix)
                expt.tix=tix;
            end
            
            if ~exist('doInteractive','var')||isempty(doInteractive)
                doInteractive=true;
            end
            
            xname=validatestring(xname,['t','segment',expt.fnames_periods,expt.fnames_trace]);
            yname=validatestring(yname,[expt.fnames_periods,expt.fnames_trace]);
            
            if xname=="segment"
                if ismember(yname,expt.fnames_periods)
                    featurePlotType='periods';
                else
                    featurePlotType='trace';
                end
            else
                xp=ismember(xname,['t',expt.fnames_periods]); %valid for period, otherwise trace
                yp=ismember(yname,expt.fnames_periods);
                if xp&&yp
                    featurePlotType='periods';
                elseif ~xp&&~yp
                    featurePlotType='trace';
                else
                    error('Incompatible feature types: mixed trace/period')
                end
            end
            
            %interactive figure: set up callbacks
            if doInteractive
                figID=gcf;
                figID.Name=['Features: ',expt.filename];
                figID.NumberTitle='off';
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
                figID.KeyPressFcn=@expt.commonKeypress;
                figID.UserData={'feature',featurePlotType,xname,yname};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            switch featurePlotType
                case 'periods'
                    expt.plot_features_periods(xname,yname);
                case 'trace'
                    expt.plot_features_trace(xname,yname);
            end
            
            axis tight
            YLIM=ylim();
            ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            switch xname
                    case {'t'}
                        xlabel('t')
                        for i=1:expt.nS
                            if strcmpi(class(hl(i)),'matlab.graphics.chart.primitive.Line')
                            hl(i).LineStyle='-';
                            end
                        end
                        hold on
                        for i=2:expt.nS
                            plot(expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
                        end
                        hold off
                        xlim([expt.t(1),expt.t(end)])
                        
                    case {'segment'} %TODO: do jitter here?
                        xticks(1:expt.nS)
                        xticklabels({expt.segment.name})
                        xlim([0.5,expt.nS+0.5])
                        xlabel('segment')
                        
                otherwise
                        xlabel(xname)
            end
            
            ylabel(yname)
            
%             end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % saving/displaying results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function displayResults(expt)
            %simple table view of results to cross reference with traces
            %examined by keypress
            %
            % keypress here to switch segment?
            %  TODO: include=true/false column, select traces for result
            %  export?
            
            if isempty(expt.resultsTrace)
                warning('results table is empty, nothing to do')
                return
            end
            
%             if isempty(expt.resfig)
%                 expt.resfig=figure('Name',['Results: ',expt.filename]);
%             else
%                 figure(expt.resfig);
%             end
            
            resfig=gcf; %uses current figure, or creates one if no figs
            resfig.Name=['Result Table: ',expt.filename];
            resfig.NumberTitle='off';
            

            colnames=[{'Trace'};{'Segment'};fieldnames(expt.segment(1).features_trace)];
            
            uit=uitable(resfig,'Data',expt.resultsTrace{:,:});
            uit.ColumnName=colnames;
            uit.RowName=[];
            uit.Units='normalized';
            uit.Position=[0.025,0.025,0.95,0.95];
            
        end
        
%         function displayResultsKeypress
%         end
        
        
        function writeToExcel(expt,outfilename,doDistributions)
            %write a header region with file/experiment info
            %default: trace level (means/stdevs) features
            
            %trace features: column 1=traceID, col2=segment, then features
            
            %option to instead do period-wise distributions - how to handle
            %jagged arrays? one segment per sheet?
            
            if isempty(expt.resultsTrace)
                warning('results table is empty, nothing to do')
                return
            end
            
            
            %header?
            
            %auto-filename...?
            if ~exist('outfilename','var')||isempty(outfilename)
                [path,fname,fext]=fileparts(expt.fullfile);
                outfilename=[path,filesep,fname,'_results',fext];
                
                %add numbers??
%                 it=1;
%                 while exist(outfilename,'file')
%                     outfilename=[path,fname,'_results',it,fext];
%                     it=it+1;
%                 end
            end
            
            writetable(expt.resultsTrace,outfilename);
        end
        
        function buildResultsTable(expt)
            ID=repmat((1:expt.nX)',expt.nS,1);
            SEG=[];TAB=table;
            for i=1:expt.nS
                SEG=[SEG;i*ones(expt.nX,1)];
                TAB=[TAB;struct2table(expt.segment(i).features_trace)];
            end
            expt.resultsTrace=table();
            expt.resultsTrace.ID=ID;
            expt.resultsTrace.Segment=SEG;
            expt.resultsTrace=[expt.resultsTrace,TAB];
        end
        
        
        function set.tix(expt,newval)
            expt.tix=newval;
            expt.updatePlots();
        end
        
    end
    
    methods (Access=private)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function detect_periods(expt,method,varargin)
            
            extras={};
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                tt=expt.t(thisIx);
                xx=expt.Xfilt(thisIx,:);
                switch method
                    case {'peaks'}
                        delta=varargin{1};
                        if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=peak_detector(tt,xx,delta,extras{:});

                    case {'threshold'}
                        frac=varargin{1};
                        if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=plateau_detector(tt,xx,frac,extras{:});
                end
                
                expt.segment(i).points=points;
                expt.segment(i).features_periods=Fdist;
                expt.segment(i).features_trace=F;
            end
            expt.featureFcn=fcns.compute_features;
        end
        
        function periodogram(expt)
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                [PSD,F,Pmax,fmax]=powerSpectrum(expt.Xdetrend(thisIx,:),expt.fs);
%                 [PSD,F,Pmax,fmax]=powerSpectrum(expt.Xfilt(thisIx,:),expt.fs);
                Tpsd=1./fmax;
                expt.psd(:,:,i)=PSD;
                
                Tpsd=num2cell(Tpsd);
                fmax=num2cell(fmax);
                Pmax=num2cell(Pmax);
                [expt.segment(i).features_trace.Tpsd]=Tpsd{:};
                [expt.segment(i).features_trace.fmax]=fmax{:};
                [expt.segment(i).features_trace.Pmax]=Pmax{:};
            end
            expt.f=F; %frquency vector
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_t(expt,whichPlot,showPts)
            
            %BUG: matlab errors when mousing over data. datatips?? super
            %annoying
            
            cla %do we want to clear current axis?  should actually behave like normal plot - add to a fig, unless undesired.
            hold on
            
            points_option='period';
            x2=[];
            switch whichPlot
                case {'r','raw'}
                    x=expt.X(:,expt.tix);
                    
                case {'n','norm','normalize','normalized'}
                    x=expt.Xnorm(:,expt.tix);
%                     x2=expt.Xtrend(:,expt.tix);

                case {'t','trend','trendline'}
                    x=expt.Xtrend(:,expt.tix);
                    
                case {'d','detrend','detrended'}
                    x=expt.Xdetrend(:,expt.tix);
%                     x2=expt.Xfilt(:,expt.tix);
                    
                case {'f','filt','filter','filtered'}
                    x=expt.Xfilt(:,expt.tix);
                    points_option='all';
                    
                otherwise
                    error([whichPlot, ' is not a supported trace to plot']);
            end
            
            plot(expt.t,x,'k')
            if ~isempty(x2)
                plot(expt.t,x2)
            end
            
            
            %slower, but benefit of encapsulating special point plots in
            %detector code
%             plotTrace=false;
%             for i=1:expt.nS
%                 tt=expt.t(expt.segment(i).ix);
%                 xx=x(expt.segment(i).ix);
%                 pts=expt.segment(i).points(expt.tix);
%                 pts.period.x=interp1(tt,xx,pts.period.t); %interp x value of period markers
%                 
%                 expt.segment(i).plot(tt,xx,pts,1,plotTrace,points_option);
%             end

            if showPts
                for i=1:expt.nS
                    tt=expt.t(expt.segment(i).ix);
                    xx=x(expt.segment(i).ix);
                    pts=expt.segment(i).points(expt.tix);
                    
                    if length(pts.period.t)>1
                        
                    if points_option=="all"
                        plot(pts.min.t,pts.min.x,'r^')
                        plot(pts.max.t,pts.max.x,'rv')
                        plot(pts.dxmin.t,pts.dxmin.x,'g<')
                        plot(pts.dxmax.t,pts.dxmax.x,'g>')
                        plot(pts.up.t,pts.up.x,'bd')
                        plot(pts.down.t,pts.down.x,'bo')
%                         line(pts.min.t,pts.min.x,'color','r','marker','^','linestyle','none')
%                         line(pts.max.t,pts.max.x,'color','r','marker','v','linestyle','none')
%                         line(pts.dxmin.t,pts.dxmin.x,'color','g','marker','<','linestyle','none')
%                         line(pts.dxmax.t,pts.dxmax.x,'color','g','marker','>','linestyle','none')
%                         line(pts.up.t,pts.up.x,'color','b','marker','d','linestyle','none')
%                         line(pts.down.t,pts.down.x,'color','b','marker','o','linestyle','none')
                        
%                         tupdwn=[pts.up.t(1:length(pts.down.t)); pts.down.t];
%                         ythresh=[1;1]*expt.segment(i).features(expt.tix).pthresh;
%                         plot(ax,tupdwn,ythresh,'b-')
                    else
                    
                        xPer=interp1(tt,xx,pts.period.t); %this is necessary except for filt
                        plot(pts.period.t,xPer,'bs')
%                         line(pts.period.t,xPer,'color','b','marker','s','linestyle','none')

                    end
                    
                    end
                end
            end
            
%             axis tight
            xlabel('t')
            ylabel(['X_{',whichPlot,'}',num2str(expt.tix)])
%             YLIM=ylim();
%             ylim([YLIM(1)-0.05*abs(YLIM(1)),YLIM(2)+0.05*abs(YLIM(2))]);
            
            for i=2:expt.nS
                plot(expt.segment(i).endpoints(1)*[1,1],ylim(),'g')
            end
            
            hold off
        end
        
        function plot_psd(expt)
            
%             axes(ax)
            cla
            hold on;
            
            for i=1:expt.nS
                plot(expt.f,expt.psd(:,expt.tix,i)); 
%                 plot(expt.f,pow2db(expt.psd(:,expt.tix,i)));
            end
            % set(gca,'yscale','log')
%             fhi=find(
%             title('One Sided Power Spectral Density');       
            xlabel('frequency')         
            ylabel('power');
            axis tight
%             xlim(ax,[0,0.5])
            if expt.nS>1
                legend({expt.segment.name})
            end
            hold off
        end
        
        function hl=plot_features_periods(expt,xname,yname)
            
            cla
            hold on;
            
            xJitter=0.05; %make this a property of the class?
           
            XX=cell(1,expt.nS);
            YY=cell(1,expt.nS);
            colorBySegment=false;  %TODO: always color by segment?
            for i=1:expt.nS
                        
                y=expt.segment(i).features_periods(expt.tix).(yname);
                if isempty(y)
                    return
                end
                
                switch xname
                    case {'t'}
                        %TODO: option for marker location?
                         x=(expt.segment(i).points(expt.tix).period.t(1:end-1)+expt.segment(i).points(expt.tix).period.t(2:end))/2; %midpoint of period
%                         x=expt.segment(i).points(expt.tix).max.t;
%                         x=expt.segment(i).points(expt.tix).min.t(1:end-1); %first min
%                         x=expt.segment(i).points(expt.tix).min.t(2:end); %second min
%                         x=expt.segment(i).points(expt.tix).down.t;
%                         x=expt.segment(i).points(expt.tix).up.t;
%                         x=(expt.segment(i).points(expt.tix).down.t+expt.segment(i).points(expt.tix).period.t(2:end))/2; %midpoint of silent post active
%                         x=(expt.segment(i).points(expt.tix).up.t+expt.segment(i).points(expt.tix).period.t(1:end-1))/2; %midpoint of silent pre active   
                        
                    case {'segment'}
                        %TODO: option for violin plot, boxplot?
                        x=i+xJitter*randn(1,length(y));

                    case expt.fnames
                        x=expt.segment(i).features(expt.tix).(xname);
                        colorBySegment=true;

                    otherwise
                        error(['Invalid name for x-axis: ' xname])
                end
                
                if isempty(x)
                    y=[];
                end
                
                XX{i}=x;
                YY{i}=y;

            end
                          
            
            ax=gca;
            ax.ColorOrderIndex=1;  
            hl=matlab.graphics.chart.primitive.Line.empty(expt.nS,0);
            for i=1:expt.nS
                if ~isempty(XX{i})
                    hl(i,1)=plot(XX{i},YY{i},'o');
                end
            end
            if colorBySegment
                legend(hl,{expt.segment.name}) %TODO: option to suppress legend?
            else
                for i=1:expt.nS
                    if ~isempty(XX{i})
                        hl(i).Color='k';
                    end
                end
            end
            
            hold off
            
        end
        
        function hl=plot_features_trace(expt,xname,yname)
            
%             axes(ax)
            cla
            hold on;
            
            xJitter=0.05;
           
            HX=[];
            HY=[];
            XX=[];
            YY=[];
            
            colorBySegment=false;
            for i=1:expt.nS
                
                switch xname
                        
                    case {'segment'}
                        %plot all trace data as line seg, highlight point expt.tix
%                         xx=i+xJitter*randn(1,expt.nX);
                        xx=i*ones(1,expt.nX);
                        yy=[expt.segment(i).features_trace.(yname)];
                        HX(end+1,1)=i;
                        HY(end+1,1)=yy(expt.tix);
                        

                    case expt.fnames_trace
                        xx=[expt.segment(i).features_trace.(xname)];
                        yy=[expt.segment(i).features_trace.(yname)];
                        colorBySegment=true;
                        HX(end+1,1)=xx(expt.tix);
                        HY(end+1,1)=yy(expt.tix);

                    otherwise
                        error(['Invalid name for x-axis: ' xname])
                end
                
                XX(end+1,:)=xx;
                YY(end+1,:)=yy;

            end
            
            %line segment below markers:
            hlAll=[];
            hlTix=[];
            if expt.nS>1
                hlAll=plot(XX,YY,'-','color',[0.5,0.5,0.5]);
                hlTix=plot(HX,HY,'k-','linewidth',1.5);
            end
            
            %markers
            
            ax=gca;
            ax.ColorOrderIndex=1;
            hmAll=plot(XX',YY','o');
            for i=1:expt.nS
                hmTix(i)=plot(HX(i),HY(i),'o');
                hmTix(i).Color=hmAll(i).Color;
                hmTix(i).MarkerFaceColor=hmAll(i).Color;
            end
            if expt.nS>1
                legend(hmAll,{expt.segment.name})
            end
            
            hold off
            
            hl.hlAll=hlAll;
            hl.hlTix=hlTix;
            hl.hmAll=hmAll;
            hl.hmTix=hmTix;
        end
        
        function updatePlots(expt)
            
            for i=1:length(expt.fig_handles)
                figure(expt.fig_handles(i));
                figData=expt.fig_handles(i).UserData;
                plotType=figData{1};
                switch plotType
                    case 'trace'
                        whichPlot=figData{2};
                        showPts=figData{3};
                        
                        nPlots=length(whichPlot);
                        for i=1:nPlots
                            subplot(nPlots,1,i)
                            expt.plot_t(whichPlot{i},showPts)
                        end
%                         expt.plot_t(whichPlot,showPts)

                    case 'psd'
                        expt.plot_psd()

                    case 'feature'
                        featurePlotType=figData{2};
                        xname=figData{3};
                        yname=figData{4};
                        switch featurePlotType
                            case 'periods'
                                expt.plot_features_periods(xname,yname);
                            case 'trace'
                                expt.plot_features_trace(xname,yname);
                        end
                end
            end
            figure(expt.active_fig);
        end
        
        function commonKeypress(expt,src,event)
            expt.active_fig=src;
            switch(event.Key)
                case {'leftarrow'}
                    if expt.tix>1
                        expt.tix=expt.tix-1;
                    end
                case {'rightarrow'}
                    if expt.tix<expt.nX
                        expt.tix=expt.tix+1;
                    end
            end
        end
        
        function figureCloseFcn(expt,src,event)
            me=gcf;
            expt.fig_handles(ismember(expt.fig_handles,me))=[];
            delete(me)
        end
    
    end
end