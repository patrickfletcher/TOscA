classdef Experiment < handle
    %class to represent a single experiment containing >=1 timeseries that share the same time sample points. For
    %example, a fluorescence videomicroscopy experiment with multiple ROIs in a field of view. 

    
    %TODO: make parameterDialog a regular figure with "go" button, like
    %selectFeaturesPopup
    
    
    %TODO: add normalization to parameterDialog
    
    %TODO: get valid options for preprocessing functions exported from those
    %functions to display in dropdown lists
    
    %TODO: forced order, but flags to choose which steps (if any) to apply?
    %     - alt: arbitrary list of operations, with fully qualified options.
    %TODO: Customizeable pipeline? eg. array processingStep objects to define
    %steps (name+operation+method+params), and pointer to which one
    %gets periods measured and which one features get measured. pages
    %of X correspond to each step: eg. X(:,:,1)=raw, X(:,:,2)=norm, ..., X(:,:,end)=per
    
    
    %TODO: plotting alternatives - plot all timeseries onces, and use an
    %update function to set the visibility on/off (also would support
    %option to plot out-of-focus traces as light gray)
    
    %TODO: updatePlots needs to be axis-aware? store plottype & extra info
    %in axis userdata not figure? store axes handles instead of figure
    %handles? => single figure GUI 
    
    %TODO: plot mean+/- stdev for per-trace features, when these are distrubtions for the trace. or boxplot. or violin.
    %if so, remove the stdev from features to plot list.
    
    
    %TODO: finalize private/public properties and methods for full OOP
        
    
    %TODO: support 3D array - 3rd dim is observable id (one [nT x nX] page per observable)
    % this requires supporting different pipeline options for each
    % observable and allows computation of new higher order features,
    % such as phase relationship
    % --> actually, will likely determine period markers using only one
    % observable, then compute features relative to those markers in
    % other traces?
    
    properties
        name=''
        date=''
        sex=''
        condition='' %eg. WT vs KO
        filename=''
        fullfilename=''
        notes={}
        
        %segment is a subinterval of the trace, the basic unit within which
        %periodicity will be detected and features will be computed
        nS
        segment=struct('name','','endpoints',[],'ix',[],...
            'points',struct(),'features_periods',struct(),'features_trace',struct())
        
        t
        dt
        fs %1/dt
        f
        psd
        
        X
        nX
        includeT %place to store trace-wise include vector
        
        groupT 
        groupID
        nG
        includeG %groupwise include vector
        
        groupMode=false
        
        featureFcn %not used?
        fnames_periods
        fnames_trace
        
        resultsTrace=table
        resultsPeriod=table
        
        
        %private (?):
        N
        group
        include
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        
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
        perSegment=false
        flattenMethod='none'
        
        filterMethod='none'
        filterParam=[]
        
        averageMethod='arithmetic'
        averageParam=[]
        
        featureMethod='threshold'
        featureParam={}
        featureExtras={}
        
        interpMethod='linear'
        
        %interactive plots:
        fig_handles=matlab.ui.Figure.empty %array of figure handles spawned by this object - use to synchronize trace in focus?
        tix %in focus
        active_fig %index into fig_handles to maintain its focus upon updates
        featurePlotType='per-period'
        xfeature='segment'
        yfeature='period'
        
        resfig=matlab.ui.Figure.empty %figure to show the results table
        featureSelectDlg
        paramDlg 
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
            fullfilename=[];
            filename=[];
            if ~isempty(varargin) && ( ischar(varargin{1}) || isstring(varargin{1}) )
                [path,filename,ext]=fileparts(varargin{1});
                fullfilename=[path,filesep,filename,ext];
            elseif ~isempty(varargin) && isnumeric(varargin{1}) && ~isempty(varargin{1})
                data=varargin{1};
                loadFile=false;
            end
            
            %load a new file
            if loadFile
                if isempty(fullfilename) || ~exist(fullfilename,'file')
                    [filename,path]=uigetfile({'*.xls*'});
                    fullfilename=[path,filename];
                end
                
                [num,txt,raw]=xlsread(fullfilename);
                %extract data
%                 ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
%                 data=cell2mat(raw(ixs:end,1:size(num,2)));
%                 ixs=find(cellfun(@(x)(isnumeric(x)&&~isnan(x)),raw(:,1)),1,'first');
                data=num;
                
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
            expt.fullfilename=fullfilename;
                
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
            
            ti=0:DT:time(end); ti=ti';
            XI=zeros(length(ti),size(Xraw,2));
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
            expt.setGroup(ones(1,expt.nX));
            expt.includeT=true(1,expt.nX);
            
            expt.tix=1;
            %groupMode default is false
            expt.group=expt.groupT; 
            expt.include=expt.includeT;
            expt.N=expt.nX;
             
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
            expt.groupT=groupvar;
            expt.groupID=unique(groupvar);
            expt.nG=length(expt.groupID);
            expt.includeG=true(1,expt.nG); %could add second input to set this at same time
            
            %update the internal current group/include/N
            if expt.groupMode
                expt.group=expt.groupID;
                expt.include=expt.includeG;
                expt.N=expt.nG;
            else
                expt.group=expt.groupT;
                expt.include=expt.includeT;
                expt.N=expt.nX;
            end
        end
        
        function excludeTraces(expt,excludeIx)
            %TODO: groupMode sensitive?
            expt.includeT(excludeIx)=false;
            if expt.groupMode
                expt.include=expt.includeG;
            else
                expt.include=expt.includeT;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Preprocessing functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function setNormalizeParameters(expt,method,methodPar)
            expt.normMethod=method;
            expt.normParam=methodPar;
        end
        
        function setDetrendParameters(expt,method,methodPar,doPerSegment,flattenMethod)
            expt.trendMethod=method;
            expt.trendParam=methodPar;
            expt.perSegment=doPerSegment;
            expt.flattenMethod=flattenMethod;
        end
        
        function setFilterParameters(expt,method,methodPar)
            expt.filterMethod=method;
            expt.filterParam=methodPar;
        end
        
        
        function preprocess(expt)
            expt.normalize();
            expt.detrend();
            if expt.groupMode
                expt.averageTraces();
            end
            expt.filter();
        end
        
        
        function normalize(expt,method,methodPar)
            if nargin>1
            expt.normMethod=method;
            expt.normParam=methodPar;
            end
            expt.Xnorm=normalizeTraces(expt.t,expt.X,expt.normMethod,expt.normParam);
        end
        
        
        %detrend supports per-segment: for piecewise linear detrend.
        function detrend(expt,method,methodPar,perSegment,flattenMethod)
            
            if nargin>1
                expt.trendMethod=method;
                expt.trendParam=methodPar;
                expt.perSegment=perSegment;
                expt.flattenMethod=flattenMethod;
            end
            
            if isempty(expt.Xnorm)
                expt.normalize();
            end
            
            if expt.perSegment
                for i=1:expt.nS
                    thisIx=expt.segment(i).ix;
                    [Xdt,Xt]=detrendTraces(expt.t(thisIx),expt.Xnorm(thisIx,:),expt.trendMethod,expt.trendParam);
                    
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
                
            else
                [expt.Xdetrend,expt.Xtrend]=detrendTraces(expt.t,expt.Xnorm,expt.trendMethod,expt.trendParam);
            end
            
        end
        
        
        function averageTraces(expt) 
            %TODO: handle edge case - all traces for a group excluded 
            
            XDTg=zeros(length(expt.t),expt.nG);
            for j=1:expt.nG  
                XDTg(:,j)=mean(expt.Xdetrend(:,expt.groupT==expt.groupID(j) & expt.includeT),2);  
            end
%             emptygroup=isnan(XDTg(1,:));
%             XDTg(:,emptygroup)=[]; %remove any groups with no members
%             expt.groupID(emptygroup)=[];
            
            expt.Xdetrend=XDTg; %groupmode overwrites detrend... TODO: keep both for plotting

        end
        
        
        function filter(expt,method,methodPar)
            if nargin>1
                expt.filterMethod=method;
                expt.filterParam=methodPar;
            end
            
            if isempty(expt.Xnorm)
                expt.normalize();
            end
            if isempty(expt.Xdetrend)
                expt.detrend();
            end
            expt.Xfilt=filterTraces(expt.t,expt.Xdetrend,expt.filterMethod,expt.filterParam);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %TODO expand this - detect_periods return only points/fcns, then compute_features should apply fcn to relevant
        %traces.
        
        function setFeatureParameters(expt,method,methodpar,varargin)
            expt.featureMethod=method;
            expt.featureParam=methodpar;
            expt.featureExtras=varargin;
        end
        
        function compute_features(expt,method,methodpar,varargin)
            if nargin>1
            expt.featureMethod=method;
            expt.featureParam=methodpar;
            expt.featureExtras=varargin;
            end
            expt.detect_periods();
            expt.periodogram();
            expt.fnames_periods=fieldnames(expt.segment(1).features_periods)';
            expt.fnames_trace=fieldnames(expt.segment(1).features_trace)';
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
%                 figID.Interruptible='off';
                figID.UserData={'trace',whichPlot,showPts};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            nPlots=length(whichPlot);
            for i=1:nPlots
                subplot(nPlots,1,i)
                expt.plot_t(whichPlot{i},showPts)
            end
            expt.active_fig=gcf;
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
%                 figID.Interruptible='off';
                figID.UserData={'psd'};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
            end
            
            expt.plot_psd();
            expt.active_fig=gcf;
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
            
            %this part is a kludge
            if exist('xname','var')&&~isempty(xname)
                xname=validatestring(xname,['t','segment',expt.fnames_periods,expt.fnames_trace]);
            else
                xname=expt.xfeature;
            end
            if exist('yname','var')&&~isempty(yname)
                yname=validatestring(yname,[expt.fnames_periods,expt.fnames_trace]);
            else
                yname=expt.yfeature;
            end
            if xname=="segment"
                if ismember(yname,expt.fnames_periods)
                    expt.featurePlotType='per-period';
                else
                    expt.featurePlotType='per-trace';
                end
            else
                xp=ismember(xname,['t',expt.fnames_periods]); %valid for period, otherwise trace
                yp=ismember(yname,expt.fnames_periods);
                if xp&&yp
                    expt.featurePlotType='per-period';
                elseif ~xp&&~yp
                    expt.featurePlotType='per-trace';
                else
                    error('Incompatible feature types: mixed trace/period')
                end
            end
            
            expt.xfeature=xname;
            expt.yfeature=yname;
            
            %interactive figure: set up callbacks
            if doInteractive
                figID=gcf;
                figID.Name=['Features: ',expt.filename];
                figID.NumberTitle='off';
                figID.KeyPressFcn=@expt.commonKeypress;
%                 figID.Interruptible='off';
                figID.UserData={'feature'};
                figID.CloseRequestFcn=@expt.figureCloseFcn;
                if ~ismember(figID,expt.fig_handles)
                    expt.fig_handles(end+1)=figID; %register the new fig with expt
                end
            end
            
            switch expt.featurePlotType
                case 'per-period'
                    expt.plot_features_periods();
                case 'per-trace'
                    expt.plot_features_trace();
            end
            expt.active_fig=gcf;
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
            
            if isempty(expt.resfig)
                expt.resfig=gcf; %uses current figure, or creates one if no figs
                expt.resfig.Name=['Result Table: ',expt.filename];
                expt.resfig.NumberTitle='off';
                expt.resfig.KeyPressFcn=@expt.commonKeypress;
    %             expt.resfig.Interruptible='off';
                expt.resfig.CloseRequestFcn=@expt.resfigCloseFcn;
            else
                figure(expt.resfig)
            end
            
            colnames=[{'Trace'};{'Group'};{'Segment'};{'Include'};fieldnames(expt.segment(1).features_trace)];
%             colnames=[{'Trace'};{'Segment'};{'Include'};fieldnames(expt.segment(1).features_trace)];
%             colformat=repmat({'numeric'},1,length(colnames));
%             coledit=false(1,length(colnames));
%             colformat(3)={'logical'}; 
%             coledit(3)=true; 
            %needs a callback function to link edited values to expt.include
            
            uit=uitable(expt.resfig,'Data',expt.resultsTrace{:,:},'tag','oscar_resultsTable');
            uit.ColumnName=colnames;
%             uit.ColumnFormat=colformat;
%             uit.ColumnEditable=coledit;
            uit.RowName=[];
            uit.Units='normalized';
            uit.Position=[0.025,0.025,0.95,0.95];
            
            figure(expt.resfig) %bring focus back to figure (for keypressfcn)
            figure(expt.active_fig) %bring focus back to active figure
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
                [path,file,fext]=fileparts(expt.fullfilename);
                outfilename=fullfile(path,[file,'_results',fext]);
                
                %suggest file name and request confirmation:
                filter='*.xls*';
                title='Enter the output file name';
                [file,path]  = uiputfile(filter,title,outfilename);
                if isequal(file,0) || isequal(path,0)
                    return
                else
                    outfilename=fullfile(path,file);
                end
            end
            
            writetable(expt.resultsTrace,outfilename);
        end
        
        function buildResultsTable(expt)
            ID=repmat((1:expt.N)',expt.nS,1);
            INC=repmat(expt.include',expt.nS,1);
            GRP=repmat(expt.group',expt.nS,1);
            SEG=[];TAB=table;
            for i=1:expt.nS
                SEG=[SEG;i*ones(expt.N,1)];
                TAB=[TAB;struct2table(expt.segment(i).features_trace)];
            end
            expt.resultsTrace=table();
            expt.resultsTrace.ID=ID;
            expt.resultsTrace.Group=GRP;
            expt.resultsTrace.Segment=SEG;
            expt.resultsTrace.Include=INC;
            expt.resultsTrace=[expt.resultsTrace,TAB];
        end
        
        
%         function set.tix(expt,newval)
%             expt.tix=newval;
%             expt.updatePlots();
%         end
        
    end
    
    methods (Access=private)
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Analysis functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function detect_periods(expt)
%             expt.featureMethod=method;
%             expt.featureParam=methodpar;
%             expt.featureExtras=varargin{:};
            
            extras=expt.featureExtras;
            for i=1:expt.nS
                thisIx=expt.segment(i).ix;
                tt=expt.t(thisIx);
                xx=expt.Xfilt(thisIx,:);
                switch expt.featureMethod
                    case {'peaks'}
                        delta=expt.featureParam;
%                         if length(varargin)>1, extras=varargin(2:end); end
                        [F, Fdist, points, fcns]=peak_detector(tt,xx,delta,extras{:});

                    case {'threshold'}
                        frac=expt.featureParam;
%                         if length(varargin)>1, extras=varargin(2:end); end
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
                
                expt.psd=zeros([size(PSD),expt.nS]);
                expt.psd(:,:,i)=PSD;
                
                PSD(1,:)=Pmax;
                [~,~,points]=peak_detector(F,PSD,0.05);
                m2=zeros(1,length(points));
                for j=1:length(points)
                    othermax=points(j).max.x~=Pmax(j);
                    if any(othermax)
                        m2(j)=max(points(j).max.x(othermax));
                    end
                end
                Rp21=m2./Pmax;
                
                Tpsd=num2cell(Tpsd);
                fmax=num2cell(fmax);
                Pmax=num2cell(Pmax);
                Rp21=num2cell(Rp21);
                
                [expt.segment(i).features_trace.Tpsd]=Tpsd{:};
                [expt.segment(i).features_trace.fmax]=fmax{:};
                [expt.segment(i).features_trace.Pmax]=Pmax{:};
                [expt.segment(i).features_trace.Rp21]=Rp21{:};
            end
            expt.f=F; %frquency vector
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function plot_t(expt,whichPlot,showPts)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            points_option='period';
            x=[];x2=[];
            switch whichPlot
                case {'r','raw'}
                    if ~expt.groupMode
                        x=expt.X(:,expt.tix);
                    end
                    
                case {'n','norm','normalize','normalized'}
                    if ~expt.groupMode
                        x=expt.Xnorm(:,expt.tix);
                        x2=expt.Xtrend(:,expt.tix);
                    end

                case {'t','trend','trendline'}
                    x=expt.Xtrend(:,expt.tix);
                    
                case {'d','detrend','detrended'}
                    x=expt.Xdetrend(:,expt.tix);
                    x2=expt.Xfilt(:,expt.tix);
                    line([min(expt.t),max(expt.t)],[0,0],'color','k','linestyle','--','tag','oscar_line')
                    
                case {'f','filt','filter','filtered'}
                    x=expt.Xfilt(:,expt.tix);
                    line([min(expt.t),max(expt.t)],[0,0],'color','k','linestyle','--','tag','oscar_line')
                    points_option='all';
                    
                otherwise
                    error([whichPlot, ' is not a supported trace type for plotting']);
            end
            
            if ~isempty(x2)
                line(expt.t,x2,'tag','oscar_line')
            end
            
            if ~isempty(x)
                line(expt.t,x,'color','k','tag','oscar_line')


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
                            line(pts.min.t,pts.min.x,'color','r','marker','^','linestyle','none','tag','oscar_line')
                            line(pts.max.t,pts.max.x,'color','r','marker','v','linestyle','none','tag','oscar_line')
                            line(pts.dxmin.t,pts.dxmin.x,'color','g','marker','<','linestyle','none','tag','oscar_line')
                            line(pts.dxmax.t,pts.dxmax.x,'color','g','marker','>','linestyle','none','tag','oscar_line')
                            line(pts.up.t,pts.up.x,'color','b','marker','d','linestyle','none','tag','oscar_line')
                            line(pts.down.t,pts.down.x,'color','b','marker','o','linestyle','none','tag','oscar_line')

    %                         tupdwn=[pts.up.t(1:length(pts.down.t)); pts.down.t];
    %                         ythresh=[1;1]*expt.segment(i).features(expt.tix).pthresh;
    %                         plot(ax,tupdwn,ythresh,'b-')
                        else

                            xPer=interp1(tt,xx,pts.period.t); %this is necessary except for filt
                            line(pts.period.t,xPer,'color','b','marker','s','linestyle','none','tag','oscar_line')

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
                    line(expt.segment(i).endpoints(1)*[1,1],ylim(),'color','g','tag','oscar_line')
                end
            end
        end
        
        function plot_psd(expt)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            for i=1:expt.nS
                line(expt.f,expt.psd(:,expt.tix,i),'tag','oscar_line'); 
            end     
            xlabel('frequency')         
            ylabel('power');
            axis tight
            if expt.nS>1
                legend({expt.segment.name},'AutoUpdate','off')
            end
        end
        
        function plot_features_periods(expt)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            xJitter=0.05; %make this a property of the class?
           
            XX=cell(1,expt.nS);
            YY=cell(1,expt.nS);
            colorBySegment=false;  %TODO: always color by segment?
            for i=1:expt.nS
                        
                y=expt.segment(i).features_periods(expt.tix).(expt.yfeature);
                if isempty(y)
                    return
                end
                
                switch expt.xfeature
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

                    case expt.fnames_periods
                        x=expt.segment(i).features_periods(expt.tix).(expt.xfeature);
                        colorBySegment=true;

                    otherwise
                        error(['Invalid name for x-axis: ' expt.xfeature])
                end
                
                if isempty(x)
                    y=[];
                end
                
                XX{i}=x;
                YY{i}=y;

            end
                          
            
            ax.ColorOrderIndex=1;  
            hl=matlab.graphics.chart.primitive.Line.empty(expt.nS,0);
            for i=1:expt.nS
                if ~isempty(XX{i})
                    hl(i,1)=line(XX{i},YY{i},'marker','o','linestyle','none','tag','oscar_line');
                end
            end
            if colorBySegment
                legend(hl,{expt.segment.name},'AutoUpdate','off') %TODO: option to suppress legend?
            else
                legend('off')
                for i=1:expt.nS
                    if ~isempty(XX{i})
                        hl(i).Color='k';
                    end
                end
            end
            
            axis tight
            
            switch expt.xfeature
                    case {'t'}
                        for i=1:expt.nS
                            if strcmpi(class(hl(i)),'matlab.graphics.chart.primitive.Line')
                            hl(i).LineStyle='-';
                            end
                        end
                        hold on
                        for i=2:expt.nS
                            line(expt.segment(i).endpoints(1)*[1,1],ylim(),'color','g','tag','oscar_line')
                        end
                        hold off
                        xlim([expt.t(1),expt.t(end)])
                        
                    case {'segment'} %TODO: do jitter here?
                        xticks(1:expt.nS)
                        xticklabels({expt.segment.name})
                        xlim([0.5,expt.nS+0.5])
            end
            
            xlabel(expt.xfeature)
            ylabel(expt.yfeature)
            
        end
        
        function plot_features_trace(expt)
            
            ax=gca;
            reset(ax)
            delete(findobj(ax,'tag','oscar_line'));
            
            HX=[];
            HY=[];
            XX=[];
            YY=[];
            
            for i=1:expt.nS
                
                switch expt.xfeature
                        
                    case {'segment'}
                        %plot all trace data as line seg, highlight point expt.tix
                        xx=i*ones(1,expt.N);
                        yy=[expt.segment(i).features_trace.(expt.yfeature)];
                        HX(end+1,1)=i;
                        HY(end+1,1)=yy(expt.tix);
                        

                    case expt.fnames_trace
                        xx=[expt.segment(i).features_trace.(expt.xfeature)];
                        yy=[expt.segment(i).features_trace.(expt.yfeature)];
                        HX(end+1,1)=xx(expt.tix);
                        HY(end+1,1)=yy(expt.tix);

                    otherwise
                        error(['Invalid name for x-axis: ' expt.xfeature])
                end
                
                XX(end+1,:)=xx;
                YY(end+1,:)=yy;

            end
            
            %line segment below markers:
            if expt.nS>1
                line(XX,YY,'color',[0.5,0.5,0.5],'linestyle','-','tag','oscar_line');
                line(HX,HY,'color','k','linestyle','-','linewidth',1.5,'tag','oscar_line');
            end
            
            %markers
            
            ax.ColorOrderIndex=1;
            hmAll=line(XX',YY','marker','o','linestyle','none','tag','oscar_line');
            axis tight
            
            for i=1:expt.nS
                line(HX(i),HY(i),'marker','o','linestyle','none',...
                    'color',hmAll(i).Color,'markerfacecolor',hmAll(i).Color,'tag','oscar_line');
            end
            
            line(XX(~expt.include)',YY(~expt.include)','color','k',...
                'linestyle','none','marker','x','tag','oscar_line');
            
            if expt.nS>1
                legend(hmAll,{expt.segment.name},'AutoUpdate','off')
            end
            
            if expt.xfeature=="segment"
                xticks(1:expt.nS)
                xticklabels({expt.segment.name})
                xlim([0.5,expt.nS+0.5])
            end
            
            if (expt.xfeature=="periodMean"||expt.xfeature=="Tpsd") && (expt.yfeature=="periodMean"||expt.yfeature=="Tpsd")
                line(xlim,xlim,'tag','oscar_line')
            end
            
            xlabel(expt.xfeature)
            ylabel(expt.yfeature)
            
        end
        
        function updatePlots(expt,type)
            
            if ~exist('type','var')
                type='all';
            end
            
            if ~isempty(expt.fig_handles)
                
            for i=1:length(expt.fig_handles)
                figData=expt.fig_handles(i).UserData;
                plotType=figData{1};
                
                if isequal(plotType,type)||type=="all"
                    figure(expt.fig_handles(i));
                    switch plotType
                        case 'trace'
                            whichPlot=figData{2};
                            showPts=figData{3};

                            nPlots=length(whichPlot);
                            for j=1:nPlots
                                %TODO: use tight_subplot?
                                subplot(nPlots,1,j);
                                expt.plot_t(whichPlot{j},showPts)
                            end

                        case 'psd'
                            expt.plot_psd()

                        case 'feature'
                            
                            switch expt.featurePlotType
                                case 'per-period'
                                    expt.plot_features_periods();
                                case 'per-trace'
                                    expt.plot_features_trace();
                            end
                    end
                end
            end
            
%             if expt.active_fig~=expt.resfig
                figure(expt.active_fig);
%             else
%                 figure(expt.fig_handles(1))
%             end
            end
            
        end
        
        function commonKeypress(expt,src,event)
            expt.active_fig=src;
            switch(event.Key)
                case {'leftarrow'}
                    %previous trace
                    if expt.tix>1
                        expt.tix=expt.tix-1;
                        expt.updatePlots();
                    end
                    
                case {'rightarrow'}
                    %next trace
                    if expt.tix<expt.N
                        expt.tix=expt.tix+1;
                        expt.updatePlots();
                    end
                    
                case 'f'
                    %select features to plot
                    expt.selectFeaturesPopup();
                    
                case 'g'
                    %pre-toggle actions
                    if expt.groupMode
                        expt.includeG=expt.include; %store include list
%                         expt.ixG=expt.tix; %store last group ix
                        expt.tix=find(expt.groupT==expt.tix,1,'first'); %determine new tix: first trace from currently in-focus group
                    else
                        expt.includeT=expt.include;
%                         expt.ixT=expt.tix;
                        expt.tix=expt.groupT(expt.tix); %determine new tix: group of currently in-focus trace 
                    end
                    
                    %toggle grouped mode
                    expt.groupMode=~expt.groupMode;
                    
                    %post-toggle actions
                    if expt.groupMode
                        expt.group=expt.groupID;
                        expt.include=expt.includeG;
                        expt.N=expt.nG;
                    else
                        expt.group=expt.groupT;
                        expt.include=expt.includeT;
                        expt.N=expt.nX;
                    end
                    
                    expt.preprocess();
                    expt.compute_features();
                    expt.updatePlots()
                    expt.displayResults();
                    
                case 'i'
                    %toggle include
                    expt.include(expt.tix)=~expt.include(expt.tix);
                    
                    expt.updatePlots('feature')
                    expt.buildResultsTable();
                    expt.displayResults();
                    
                case 'p'
                    %parameter dialog
                    expt.parameterDialog();
                    %needs full pipeline helper function
                    
                case 's'
                    %save results
                    expt.writeToExcel();
                    
                case 'q'
                    %quit - close all windows
                    expt.quit();
            end
        end
        
        function figureCloseFcn(expt,~,~)
            me=gcf;
            expt.fig_handles(ismember(expt.fig_handles,me))=[];
            delete(me)
        end
        
        function resfigCloseFcn(expt,~,~)
            me=gcf;
            if isequal(expt.resfig,me) 
                expt.resfig=[];
                delete(me)
            end
        end
        
        function selectFeaturesPopup(expt)
            
            expt.featureSelectDlg=figure('Name','Select features to plot','NumberTitle','off','MenuBar','none');
            expt.featureSelectDlg.KeyPressFcn=@expt.commonKeypress;
            expt.featureSelectDlg.Position(3:4)=[300,150];
            uibg=uibuttongroup(expt.featureSelectDlg,'Position',[0.1,0.75,0.8,0.2],'SelectionChangedFcn',@changeRadio);
%             uibg=uibuttongroup(figh,'Position',[0,80,300,40],'SelectionChangedFcn',@changeRadio);
            uirb1=uicontrol(uibg,'Style','radiobutton','String','per-period');
            uirb1.Units='normalized'; uirb1.Position=[0,0,0.5,1];
            uirb2=uicontrol(uibg,'Style','radiobutton','String','per-trace');
            uirb2.Units='normalized'; uirb2.Position=[0.5,0,0.5,1];
            if ~strcmp(expt.featurePlotType,uirb1.String)
                uirb2.Value=true;
            end
            
            switch expt.featurePlotType
                case 'per-period'
                    selections={[{'t'},{'segment'},expt.fnames_periods],expt.fnames_periods};
                case 'per-trace'
                    selections={[{'t'},{'segment'},expt.fnames_trace],expt.fnames_trace};
            end
            
            uit=uitable(expt.featureSelectDlg,'Units','Normalized','Position',[0.1,0.3,0.8,0.4]);
            uit.RowName=[];
            uit.ColumnWidth={118,118};
            uit.ColumnEditable=true;
            uit.ColumnName={'x feature','y feature'};
            uit.ColumnFormat=selections;
            uit.Data={expt.xfeature,expt.yfeature};
            uit.CellEditCallback=@cellEdit;
            
            if expt.featurePlotType=="per-trace"
                uibg.SelectedObject=uirb2;
            end
            
            uicontrol('Units','Normalized','Position',[0.2,0.05,0.6,0.2],'String','Plot',...
              'Callback',@launchFeatureUpdate);
          
            function launchFeatureUpdate(~,~)
                
                for i=1:length(expt.fig_handles)
                    if strcmp(expt.fig_handles(i).UserData,'feature')
                        expt.active_fig=expt.fig_handles(i);
                    end
                end
                
                expt.updatePlots('feature')
            end
            
            function changeRadio(~,cbdata)
                switch cbdata.NewValue.String
                    case 'per-period'
                        valid_features={[{'t'},{'segment'},expt.fnames_periods],expt.fnames_periods};
                        uit.ColumnFormat=valid_features;
                        if ~ismember(uit.Data(1),valid_features{1})
                            expt.xfeature=valid_features{1}{1};
                            uit.Data(1)=valid_features{1}(1);
                        end
                        if ~ismember(uit.Data(2),valid_features{2})
                            expt.yfeature=valid_features{2}{1};
                            uit.Data(2)=valid_features{2}(1);
                        end
                            
                    case 'per-trace'
                        valid_features={[{'segment'},expt.fnames_trace],expt.fnames_trace};
                        uit.ColumnFormat=valid_features;
                        if ~ismember(uit.Data(1),valid_features{1})
                            expt.xfeature=valid_features{1}{1};
                            uit.Data(1)=valid_features{1}(1);
                        end
                        if ~ismember(uit.Data(2),valid_features{2})
                            expt.yfeature=valid_features{2}{1};
                            uit.Data(2)=valid_features{2}(1);
                        end
                end
                expt.featurePlotType=cbdata.NewValue.String;
            end
            
            function cellEdit(~,cbdata)
                if cbdata.Indices(2)==1
                    expt.xfeature=cbdata.EditData;
                elseif cbdata.Indices(2)==2
                    expt.yfeature=cbdata.EditData;
                end
            end

        end
        
       
        function parameterDialog(expt)
            prompt = {'Detrend:','Filter:','Threshold:'};
            title = 'Enter parameter values';
            dims = [1,50];
            oldtrend=expt.trendParam;
            oldfilt=expt.filterParam;
            oldfeat=expt.featureParam;
            definput = {num2str(oldtrend), num2str(oldfilt), num2str(oldfeat,'%g, %g')};
            opts.WindowStyle='normal';
            answer = inputdlg(prompt,title,dims,definput,opts);
            if isempty(answer)
                return
            else
                doUpdate=false;
                newtrend = str2double(answer{1});
                if newtrend~=oldtrend
                    expt.trendParam=newtrend;
                    doUpdate=true;
                end
                
                newfilt = str2double(answer{2});
                if newfilt~=oldfilt
                    expt.filterParam=newfilt;
                    doUpdate=true;
                end
                
                newfeat = str2num(answer{3});
                if any(newfeat~=oldfeat)
                    expt.featureParam=newfeat;
                    doUpdate=true;
                end
                
                if doUpdate==true
                    expt.preprocess();
                    expt.compute_features();
                    expt.updatePlots()
                    expt.displayResults();
                end
            end
        end
        
        function quit(expt)
            if ~isempty(expt.fig_handles)
                close(expt.fig_handles)
            end
            if ~isempty(expt.resfig)
                close(expt.resfig)
            end
            if ~isempty(expt.featureSelectDlg)
                close(expt.featureSelectDlg)
            end
        end
    end
end