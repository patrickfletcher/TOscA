classdef TOscA<handle
    %TOSCA TOscA is a graphical interface for performing
    %processing and feature extraction on timeseries data. 
    
    %TOscA - Timeseries Oscillation Analysis tool
    
    %TODO: raw/normalized data - for any feature with units of X eg. ratio
    %      vs filtered/rescaled/etc - for period/plateau detection
    % - set of normalization options, set of filtering options?
    
    %data and processing properties
    properties
        
        params
        
        t
        dt
        fs
        tIntervals
        tIntOmit
        
        Xraw
        Xnorm
        Xtrend
        Xfilt
        nTraces
        
        thrPrct=[25,75]
        thrFrac=[0.5,0.4]
        Thi=1
        Tlo=20
        
        trace
        F
    end
    
    %GUI properties
    properties
        
        hFig
        hAxRaw
        hAxNorm
        hAxFilt
        
        tix=1 %focus trace
        oldix=1
        lightgray=[0.75,0.75,0.75]
        bgLineWidth=0.5
        boldLineWidth=1.5
        
        hLraw %array of line object handles (one per trace in dataset)
        hLnorm
        hLfilt
        filtRange
        
        hMarksR %markers for special points
        hMarksN
        hMarksF
        hThrUp
        hThrDown
    end
    
    %constructor
    methods
        function app=TOscA(data,tCond)
            %assume data is either filename, or data matrix [t,X], columns
            %are observations.
            
            %TODO: also assume time is in seconds - convenient way to set this?
            
            %TODO: input checking
            %if no inputs, request file. must be able to map at least to t, Xraw
            if nargin==0
                app.loadData();
                tCond=[]; %need to be able to set this
            else
                app.t=data(:,1);
                app.Xraw=data(:,2:end);
            end
            
            app.setParams();
            
            app.dt=mode(diff(app.t));
            app.fs=1/app.dt;
            
            %set up time intervals
            %TODO: adjustable endpoints - omit first bit of each interval?
            %TODO: use omitted data for smoothing, or just mirror for edge?
            tOmit=0;
            tCond=[app.t(1),tCond(:)',app.t(end)];
            for i=1:length(tCond)-1
                app.tIntervals(i,:)=[tCond(i),tCond(i+1)];
                app.tIntOmit(i,:)=[tCond(i)+tOmit,tCond(i+1)];
            end
            
            app.nTraces=size(app.Xraw,2);
            
            app.buildGUI(); %set up display elements
            app.processData(); %to get Xfilt
            app.getFeatures();
            app.plotData(); %populate the axes
        end
    end
    
    %data processing and feature detection functions
    methods
        
        %apply signal processing steps to raw data
        function processData(app)
            % TEST:
            %bandpass whole trace, amp-normalize+feature detect in IOIs?
            
            app.Xfilt=nan(size(app.Xraw));
            
            for i=1:size(app.tIntervals,1)
                ix=app.t>=app.tIntervals(i,1)&app.t<=app.tIntervals(i,2); %both have = to catch last point...
                XRAW=app.Xraw(ix,:);
                
                %normailze
                XNORM=XRAW./mean(XRAW,1);
                
                %detrend
                pars.trend.method='lowpass'; %for better filter design options than moving average
                pars.trend.fpass=1/app.Tlo;%cutoff period for lowpass to create smoothed traces as trend
                [XDT,XT]=detrend(app.t,XNORM,pars);
                
                %filter
                %lowpass only
                XFILT=lowpass(XDT,1/app.Thi,app.fs);
                
                %bandpass method 1
                %                 wszTrend=13;
                %                 wszLP=1;
                %                 wszTrend=round(wszTrend/app.dt);
                %                 wszTrend=max(wszTrend,3);
                %                 wszLP=round(wszLP/app.dt);
                %                 wszLP=max(wszLP,3);
                %                 XTrend=smoothdata(XRAW,method,wszTrend);
                %                 XHP=XRAW-XTrend;
                %                 XFILT=smoothdata(XHP,method,wszLP);
                
                %bandpass method 2
%                 T1=12; T2=1;
%                 fpass=[1/T1,1/T2]; %Hz
%                 XFILT=bandpass(XRAW,fpass,app.fs,'Steepness',0.5,'ImpulseResponse','iir');
                
                
                %final normalization useful for getting ap/sp/period
%                 maxX=max(XLP,[],1);
%                 minX=min(XLP,[],1);
%                 XFILT=(XFILT-minX)./abs(maxX-minX);

%                 p=prctile(XFILT,[25,75],1);
%                 p=prctile(XFILT,[5,95],1);
%                 p=prctile(XFILT,[2.5,97.5],1);
                p=prctile(XFILT,[0,100],1);
%                 p=prctile(XFILT,app.thrPrct,1);
                XFILT=(XFILT-p(1,:))./abs(p(2,:)-p(1,:));
                app.filtRange=p;
                
%                 meanX=mean(XFILT,1);
%                 XFILT=(XFILT-meanX);
                
                %envelope detection and amplitude normalization
%                 wszTrend=20;
%                 wszTrend=round(wszTrend/app.dt);
%                 wszTrend=max(wszTrend,3);
%                 [envMax,envMin]=envelope(XFILT,wszTrend,'rms');
%                 XFILT=(XFILT-envMin)./abs(envMax-envMin);
%                 
                
                app.Xnorm(ix,:)=XDT;
                app.Xfilt(ix,:)=XFILT;
                
                %                 ixOmit=find(app.t>=app.tIntOmit(i,1)&app.t<=app.tIntOmit(i,2));
                %                 app.Xfilt(ixOmit,:)=XLP(ixOmit,:);
            end
        end
        
        %threshold-based feature detection
        function getFeatures(app)
            %use Xfilt to find time features (up/down transitions)
            %use Xraw then for max/min in up/down phases.
            %alt: use Xfilt max/min times, but use Xraw values.
            
            %computed features per period:
            % Period=tup(i+1)-tup(i)
            % AP=tdown(i)-tup(i)
            % PlatFrac=AP/Period
            %
            % Amplitude=yMax(i)-yMin(i)
            %
            % Alt Period: from PSD
            
            %statistics:
            % - mean over total periods observed.
            % - variance propagation through computed values?
            
            %run plateau detector on full processed trace, partition events into intervals
            [Pts,F]=plateau_detector(app.t, app.Xfilt, app.thrFrac,'ThresholdPercentiles',app.thrPrct);
            for i=1:app.nTraces
                if ~isempty(Pts(i).tUp)
                Pts(i).xUpR=interp1(app.t,app.Xraw(:,i),Pts(i).tUp);
                Pts(i).xUpN=interp1(app.t,app.Xnorm(:,i),Pts(i).tUp);
                end
                if ~isempty(Pts(i).tDown)
                Pts(i).xDownR=interp1(app.t,app.Xraw(:,i),Pts(i).tDown);
                Pts(i).xDownN=interp1(app.t,app.Xnorm(:,i),Pts(i).tDown);
                end
            end
            for interval=1:size(app.tIntervals,1)
                
                   
                thisPts=Pts;
                for i=1:app.nTraces
                    ix=Pts(i).tUp>=app.tIntervals(interval,1)&Pts(i).tUp<=app.tIntervals(interval,2); %both have = to catch last point...
                    thisPts(i).tUp=Pts(i).tUp(ix);
                    thisPts(i).xUp=Pts(i).xUp(ix);
                    thisPts(i).xUpR=Pts(i).xUpR(ix);
                    thisPts(i).xUpN=Pts(i).xUpN(ix);
                    
                    ix=Pts(i).tDown>=app.tIntervals(interval,1)&Pts(i).tDown<=app.tIntervals(interval,2); %both have = to catch last point...
                    thisPts(i).tDown=Pts(i).tDown(ix);
                    thisPts(i).xDown=Pts(i).xDown(ix);
                    thisPts(i).xDownR=Pts(i).xDownR(ix);
                    thisPts(i).xDownN=Pts(i).xDownN(ix);
                end
                
                app.trace{interval}.points=thisPts;
                
%                 app.trace{interval}.features=thisF;
            end

%             %run plateau detector separately on each interval.
%             for interval=1:size(app.tIntervals,1)
%                 ix=app.t>=app.tIntervals(interval,1)&app.t<=app.tIntervals(interval,2); %both have = to catch last point...
%                 [Pts,F]=plateau_detector(app.t(ix), app.Xfilt(ix,:), [0.4,0.3]);
%                 
%                 % interpolate the values of Xraw/Xnorm at feature detector times
%                 for i=1:app.nTraces
%                     if ~isempty(Pts(i).tUp)
%                     Pts(i).xUpR=interp1(app.t(ix),app.Xraw(ix,i),Pts(i).tUp);
%                     Pts(i).xUpN=interp1(app.t(ix),app.Xnorm(ix,i),Pts(i).tUp);
%                     end
%                     if ~isempty(Pts(i).tDown)
%                     Pts(i).xDownR=interp1(app.t(ix),app.Xraw(ix,i),Pts(i).tDown);
%                     Pts(i).xDownN=interp1(app.t(ix),app.Xnorm(ix,i),Pts(i).tDown);
%                     end
%                 end
%                 app.trace{interval}.points=Pts;
%                 app.trace{interval}.features=F;
%             end
            
        end
    end
    
    %callbacks and other GUI functions
    methods
        %modify default analysis parameters
        %TODO: UI table? pops up if needed?
        function setParams(app)
        end
        
        %rerun analysis and plot data (eg after changing params)
        function updateAll(app)
            app.processData(); %to get Xfilt
            app.plotData(); %populate the axes
        end
        
        %set up all the graphical components of the UI
        function buildGUI(app)
            app.hFig=figure('Name','Trace Analyzer','NumberTitle','off');
            app.hFig.KeyPressFcn=@app.keyPressDecoder;
            
            app.hAxRaw=subplot(3,1,1); 
            app.hAxNorm=subplot(3,1,2);
            app.hAxFilt=subplot(3,1,3);
            
        end
        
        function plotData(app)
            app.hLraw=plot(app.hAxRaw,app.t,app.Xraw,'color',app.lightgray,'LineWidth',app.bgLineWidth); 
            app.hLnorm=plot(app.hAxNorm,app.t,app.Xnorm,'color',app.lightgray,'LineWidth',app.bgLineWidth);
            app.hLfilt=plot(app.hAxFilt,app.t,app.Xfilt,'color',app.lightgray,'LineWidth',app.bgLineWidth); 
            
            ylabel(app.hAxRaw,'Xraw')
            axis(app.hAxRaw,'tight')
            
            ylabel(app.hAxNorm,'Xnorm')
            axis(app.hAxNorm,'tight')
            
            xlabel(app.hAxFilt,'Time')
            ylabel(app.hAxFilt,'Xfilt')
            axis(app.hAxFilt,'tight')
            
            linkaxes([app.hAxRaw,app.hAxNorm,app.hAxFilt],'x')
            
            hold(app.hAxRaw,'on')
            hold(app.hAxNorm,'on')
            hold(app.hAxFilt,'on')
            
            app.updateCurrentTrace();
        end
        
        function updateCurrentTrace(app)
            %raw traces
            app.hLraw(app.oldix).Color=app.lightgray;
            app.hLraw(app.oldix).LineWidth=app.bgLineWidth;
            app.hLraw(app.oldix).ZData=zeros(size(app.t));
            
            app.hLraw(app.tix).Color='k'; %black
            app.hLraw(app.tix).LineWidth=app.boldLineWidth; %bold
            app.hLraw(app.tix).ZData=ones(size(app.t)); %on top of others
            
            %raw traces
            app.hLnorm(app.oldix).Color=app.lightgray;
            app.hLnorm(app.oldix).LineWidth=app.bgLineWidth;
            app.hLnorm(app.oldix).ZData=zeros(size(app.t));
            
            app.hLnorm(app.tix).Color='k'; %black
            app.hLnorm(app.tix).LineWidth=app.boldLineWidth; %bold
            app.hLnorm(app.tix).ZData=ones(size(app.t)); %on top of others
            
            %filtered traces
            app.hLfilt(app.oldix).Color=app.lightgray;
            app.hLfilt(app.oldix).LineWidth=app.bgLineWidth;
            app.hLfilt(app.oldix).ZData=zeros(size(app.t));
            
            app.hLfilt(app.tix).Color='k';
            app.hLfilt(app.tix).LineWidth=app.boldLineWidth;
            app.hLfilt(app.tix).ZData=ones(size(app.t));
            
            delete(app.hMarksR)
            delete(app.hMarksN)
            delete(app.hMarksF)
            delete(app.hThrUp)
            delete(app.hThrDown)
            
            hup=[]; hdwn=[];
            
            for interval=1:size(app.tIntervals,1)
                Pts=app.trace{interval}.points(app.tix);
                
                if ~isempty(Pts.tUp)
                    hup=plot(app.hAxRaw,Pts.tUp,Pts.xUpR,'bs');
                    app.hMarksR=[app.hMarksR,hup];
                    
                    hup=plot(app.hAxNorm,Pts.tUp,Pts.xUpN,'bs');
                    app.hMarksN=[app.hMarksN,hup];
                    
                    hup=plot(app.hAxFilt,Pts.tUp,Pts.xUp,'bs');
                    app.hMarksF=[app.hMarksF,hup];
                    thrUp=Pts.xUp(1);
                    app.hThrUp=plot(app.hAxFilt,xlim(),thrUp*[1,1],'r--');
                end
                if ~isempty(Pts.tDown)
                    hdwn=plot(app.hAxRaw,Pts.tDown,Pts.xDownR,'bo');
                    app.hMarksR=[app.hMarksR,hdwn];
                    
                    hdwn=plot(app.hAxNorm,Pts.tDown,Pts.xDownN,'bo');
                    app.hMarksN=[app.hMarksN,hdwn];
                    
                    hdwn=plot(app.hAxFilt,Pts.tDown,Pts.xDown,'bo');
                    app.hMarksF=[app.hMarksF,hdwn];
                    thrDown=Pts.xDown(1);
                    if thrUp~=thrDown
                        app.hThrDown=plot(xlim(),thrDown*[1,1],'b--');
                    end
                end

            end
            
%             ylim(app.hAxNorm,app.filtRange(:,app.tix));
%             ylim(app.hAxFilt,prctile(app.Xfilt(:,app.tix),app.thrPrct));
            
        end
        
        function keyPressDecoder(app, src, event)
            %src - object in focus when keypress occurred
            %event - details of the keypress event
            
            % uicontrol-dependent keypress hierarchy of actions
            
            switch(event.Key)
                case {'leftarrow'}
                    if app.tix>1
                        app.oldix=app.tix;
                        app.tix=app.tix-1;
                        app.updateCurrentTrace();
                    end
                case {'rightarrow'}
                    if app.tix<app.nTraces
                        app.oldix=app.tix;
                        app.tix=app.tix+1;
                        app.updateCurrentTrace();
                    end
                case {'p'}
                    disp('change params')
                case {'s'}
                    disp('save placeholder')
                case {'l'}
                    disp('load placeholder')
                    
            end
        end
        
    end
    
    %File I/O
    methods
        
        %parse input arguments
        function parseArguments(app,varargin)
            %inputparser? may not be necessary at all.
        end
        
        %select and load an excel file to read data from
        function loadData(app,filename)
            %allow header to contain parameters and metadata for the
            %experiment, eg. condition times, condition names, etc.
            
            if exist('filename','var')&&~isempty(filename)
                %try to load, if fails uigetfile
            else
                [filename,path]=uigetfile({'*.xls*'});
                filename=[path,filename];
            end
            
            %load data
            [data,header]=xlsread(filename);
            time=data(:,1);
            X=data(:,2:end);
            
            %decode header info
            
            %interpolate any missing data using neighboring time points
            rowIsNan=any(isnan(X),2);
            if any(rowIsNan)
%                 flag=[flag,{'some rows have NaN'}];
                r=find(rowIsNan);
                if r(1)==1
                    time=time(2:end,1);
                    X=X(2:end,:);
                    r=r(2:end)-1;
                end
                for ii=1:length(r)
                    for jj=1:size(X,2)
                        iix=[r(ii)-1,r(ii)+1];
                        x=time(iix,1);
                        v=X(iix,jj);
                        xq=time(r(ii),1);
                        X(r(ii),jj)=interp1(x,v,xq);
                    end
                end
            end
            
            app.t=time;
            app.Xraw=X;
        end
        
        %save results/params: MAT + XLS? ask for save name
        function saveResults(app)
            %collect results for writing to file
            %uiputfile
        end
    end
    
end