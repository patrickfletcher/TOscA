classdef OscillationAnalyzer<handle
    %OscillationAnalyzer is a graphical interface for performing
    %processing and feature extraction on timeseries data.
    
    %TOscA - Timeseries Oscillation Analysis tool
    
    
    %data and processing properties
    properties
        
        params
        
        t
        dt
        fs
        
        condTimes
        condNames
        tIntervals
        nIntervals
        tOmit=0
        nOmit=0
        
        Xraw
        Xnorm
        Xtrend
        Xdetrend
        Xfilt
        XfiltUnit
        nTraces
        
        thrPrct=[0,100]
        thrFrac=[0.5,0.45]
        Thi=1
        Tlo=20
        
        trace
    end
    
    %GUI properties
    properties
        
        datafilepath
        datafilename
        
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
        
    end
    
    %constructor
    methods
        function app=OscillationAnalyzer(data,condTimes,condNames)
            %assume data is either filename, or data matrix [t,X], columns
            %are observations.
            
            %TODO: also assume time is in seconds - convenient way to set this?
            
            %TODO: input checking
            %if no inputs, request file. must be able to map at least to t, Xraw
            if nargin==0
                app.loadData();
                condTimes=[]; %need to be able to set this
            else
                if isnumeric(data)
                    app.t=data(:,1);
                    app.Xraw=data(:,2:end);
                else
                    filename=data;
                    app.loadData(filename);
                end
                
                if ~exist('condTimes','var')
                    condTimes=[];
                end
                if ~exist('condNames','var')
                    condNames={''};
                end
                
                app.dt=mode(diff(app.t));
                app.fs=1/app.dt;
                
                %set up time intervals
                %TODO: adjustable endpoints - omit first bit of each interval?
                %TODO: use omitted data for smoothing, or just mirror for edge?
                %             tOmit=0;
                app.nOmit=round(app.tOmit/app.dt);
                condTimes=unique([app.t(1),condTimes(:)',app.t(end)]);
                app.nIntervals=length(condTimes)-1;
                for i=1:app.nIntervals
                    app.tIntervals(i,:)=[condTimes(i),condTimes(i+1)];
                end
                app.condTimes=condTimes(1:end-1);
                app.condNames=condNames;
                
                app.nTraces=size(app.Xraw,2);
            end
            
            %             app.setParams();
            app.buildGUI(); %set up display elements
            app.processData();
            app.plotData(); %populate the axes
        end
    end
    
    %data processing and feature detection functions
    methods
        
        %apply signal processing steps to raw data
        function processData(app)
            % TEST:
            %bandpass whole trace, amp-normalize+feature detect in IOIs?
            
            app.Xnorm=nan(size(app.Xraw));
            app.Xfilt=nan(size(app.Xraw));
            
            %             for i=1:size(app.tIntervals,1)
            %                 ix=app.t>=app.tIntervals(i,1)&app.t<=app.tIntervals(i,2); %both have = to catch last point...
            
            ix=true(size(app.t));
            %
            T=app.t(ix);
            XRAW=app.Xraw(ix,:);
            
            %normailze
%                             method='none';
%             method='devmean';
            method='devmean2';
            %                 method='devmedian';
            %                 method='zscore';
            %                 method='unit';
            methodpar=[];
            
            %                 method='scale';
            %                 methodpar='mean';
            % %                 methodpar='median';
            
            %                 method='trend';
            % %                 methodpar={'gaussian',app.Tlo};
            %                 methodpar={'lowpass',1/app.Tlo};
            XNORM=normalizeTraces(T,XRAW,method,methodpar);
            
            %detrend
%                             method='none';
%                             methodpar=[];
            %                 method='gaussian';
            method='sgolay';
            methodpar=app.Tlo;
            %                 method='lowpass';
            %                 methodpar=1/app.Tlo;%cutoff period for lowpass to create smoothed traces as trend
            [XDT,XT]=detrendTraces(T,XNORM,method,methodpar);
            
            %filter
            %                 method='none';
            %                 methodpar=[];
                            method='gaussian'; %for better filter design options than moving average
                            methodpar=app.Thi;
%             method='lowpass'; %for better filter design options than moving average
%             methodpar=1/app.Thi;
            XFILT=filterTraces(T,XDT,method,methodpar);
            
            
            %final normalization per interval, useful for threshold detect
            XFILT=normalizeTraces(T,XFILT,'wptile',{15,[0,100]});
            XFILT(XFILT>1)=1;
            XFILT(XFILT<0)=0;

%             for i=1:size(app.tIntervals,1)
%                 ix=app.t>=app.tIntervals(i,1)&app.t<=app.tIntervals(i,2); %both have = to catch last point...
% %                 XFILT(ix,:)=normalizeTraces(T(ix),XFILT(ix,:),'unit',[]);
%                 XFILT(ix,:)=normalizeTraces(T(ix),XFILT(ix,:),'ptile',[5,95]);
% 
%             end
            
            %envelope detection and amplitude normalization
            %                 wszTrend=20;
            %                 wszTrend=round(wszTrend/app.dt);
            %                 wszTrend=max(wszTrend,3);
            %                 [envMax,envMin]=envelope(XFILT,wszTrend,'rms');
            %                 XFILT=(XFILT-envMin)./abs(envMax-envMin);
            %
            
            %                 app.Xnorm(ix,:)=XDT;
            %                 app.Xfilt(ix,:)=XFILT;
            %
            %                 [P,f,Pmax,fmax]=powerSpectrum(XFILT-mean(XFILT,1),app.fs);
            %                 app.trace{i}.features.Tpsd=1./fmax;
            %                 app.trace{i}.features.Ppsd=Pmax;
            %             end
            
            app.Xnorm=XNORM;
            app.Xtrend=XT;
            app.Xdetrend=XDT;
            app.Xfilt=XFILT;
            
            app.getFeatures();
            
        end
        
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
            %             [Pts,F]=plateau_detector(app.t, app.Xfilt, app.thrFrac);
            % %             [Pts,F]=plateau_detector(app.t, app.Xfilt, app.thrFrac,'ThresholdPercentiles',app.thrPrct);
            %             for i=1:app.nTraces
            %                 if ~isempty(Pts(i).tUp)
            %                 Pts(i).xUpR=interp1(app.t,app.Xraw(:,i),Pts(i).tUp);
            %                 Pts(i).xUpN=interp1(app.t,app.Xnorm(:,i),Pts(i).tUp);
            %                 end
            %                 if ~isempty(Pts(i).tDown)
            %                 Pts(i).xDownR=interp1(app.t,app.Xraw(:,i),Pts(i).tDown);
            %                 Pts(i).xDownN=interp1(app.t,app.Xnorm(:,i),Pts(i).tDown);
            %                 end
            %             end
            %             for interval=1:size(app.tIntervals,1)
            %
            %
            %                 thisPts=Pts;
            %                 for i=1:app.nTraces
            %                     ix=Pts(i).tUp>=app.tIntervals(interval,1)&Pts(i).tUp<=app.tIntervals(interval,2); %both have = to catch last point...
            %                     thisPts(i).tUp=Pts(i).tUp(ix);
            %                     thisPts(i).xUp=Pts(i).xUp(ix);
            %                     thisPts(i).xUpR=Pts(i).xUpR(ix);
            %                     thisPts(i).xUpN=Pts(i).xUpN(ix);
            %
            %                     ix=Pts(i).tDown>=app.tIntervals(interval,1)&Pts(i).tDown<=app.tIntervals(interval,2); %both have = to catch last point...
            %                     thisPts(i).tDown=Pts(i).tDown(ix);
            %                     thisPts(i).xDown=Pts(i).xDown(ix);
            %                     thisPts(i).xDownR=Pts(i).xDownR(ix);
            %                     thisPts(i).xDownN=Pts(i).xDownN(ix);
            %                 end
            %
            %                 app.trace{interval}.points=thisPts;
            %
            % %                 app.trace{interval}.features=thisF;
            %             end
            
            %run plateau detector separately on each interval.
            for interval=1:size(app.tIntervals,1)
                ix=app.t>=app.tIntervals(interval,1)&app.t<=app.tIntervals(interval,2); %both have = to catch last point...
                T=app.t(ix);
                XRAW=app.Xraw(ix,:);
                XNORM=app.Xnorm(ix,:);
                XFILT=app.Xfilt(ix,:);
                [Pts,F]=plateau_detector(T, XFILT, app.thrFrac);
                
                % interpolate the values of Xraw/Xnorm at feature detector times
                for i=1:app.nTraces
                    nT=length(Pts(i).tUp)-1;
                    if nT>0
                        Pts(i).xUpR=interp1(T,XRAW(:,i),Pts(i).tUp);
                        Pts(i).xUpN=interp1(T,XNORM(:,i),Pts(i).tUp);
                        
                        Pts(i).xDownR=interp1(T,XRAW(:,i),Pts(i).tDown);
                        Pts(i).xDownN=interp1(T,XNORM(:,i),Pts(i).tDown);
                        
                        for j=1:nT
                            tt=Pts(i).iUp(j):Pts(i).iUp(j+1)-1;
                            [xmax,imax]=max(XRAW(tt,i));
                            Pts(i).iMaxR(j)=tt(1)+imax-1;
                            Pts(i).tMaxR(j)=T(tt(imax));
                            Pts(i).xMaxR(j)=xmax;
                            
                            [xmax,imax]=max(XNORM(tt,i));
                            Pts(i).iMaxN(j)=tt(1)+imax-1;
                            Pts(i).tMaxN(j)=T(tt(imax));
                            Pts(i).xMaxN(j)=xmax;
                            
                            [xmin,imin]=min(XRAW(tt,i));
                            Pts(i).iMinR(j)=tt(1)+imin-1;
                            Pts(i).tMinR(j)=T(tt(imin));
                            Pts(i).xMinR(j)=xmin;
                            
                            [xmin,imin]=min(XNORM(tt,i));
                            Pts(i).iMinN(j)=tt(1)+imin-1;
                            Pts(i).tMinN(j)=T(tt(imin));
                            Pts(i).xMinN(j)=xmin;
                        end
                    end
                end
                app.trace{interval}.points=Pts;
                app.trace{interval}.features=F;
            end
            
        end
    end
    
    %callbacks and other GUI functions
    methods
        %modify default analysis parameters
        %TODO: UI table? pops up if needed?
        function setParams(app)
        end
        
        function updateTimeIntervals(app)
        end
        
        %set up all the graphical components of the UI
        function buildGUI(app)
            app.hFig=figure('Name',['Oscillation Analyzer - ',app.datafilename],'NumberTitle','off');
            app.hFig.KeyPressFcn=@app.keyPressDecoder;
%             app.hFig.Position=[20,20,];
            
            app.hAxRaw=subplot(3,1,1);
            app.hAxNorm=subplot(3,1,2);
            app.hAxFilt=subplot(3,1,3);
            
        end
        
        function plotData(app)
            
            app.hLraw=plot(app.hAxRaw,app.t,app.Xraw,'color',app.lightgray,'LineWidth',app.bgLineWidth);
            app.hLnorm=plot(app.hAxNorm,app.t,app.Xnorm,'color',app.lightgray,'LineWidth',app.bgLineWidth);
            app.hLfilt=plot(app.hAxFilt,app.t,app.Xfilt,'color',app.lightgray,'LineWidth',app.bgLineWidth);
            
            Tlim=[min(app.t),max(app.t)];
            
            ylabel(app.hAxRaw,'Xraw')
            %             axis(app.hAxRaw,'tight')
            app.hAxRaw.XLim=Tlim;
            app.hAxRaw.YLim=[min(app.Xraw(:)),max(app.Xraw(:))];
            
            ylabel(app.hAxNorm,'Xnorm,detrend')
            %             axis(app.hAxNorm,'tight')
            app.hAxNorm.XLim=Tlim;
            app.hAxNorm.YLim=[min(app.Xnorm(:)),max(app.Xnorm(:))];
            
            xlabel(app.hAxFilt,'Time')
            ylabel(app.hAxFilt,'Xfilt')
            %             axis(app.hAxFilt,'tight')
            app.hAxFilt.XLim=Tlim;
            app.hAxFilt.YLim=[min(app.Xfilt(:)),max(app.Xfilt(:))];
            
            for i=1:app.nIntervals-1
                line(app.hAxRaw,app.tIntervals(i,2)*[1,1],app.hAxRaw.YLim,'color',[0,0.75,0]);
                line(app.hAxNorm,app.tIntervals(i,2)*[1,1],app.hAxNorm.YLim,'color',[0,0.75,0]);
                line(app.hAxFilt,app.tIntervals(i,2)*[1,1],app.hAxFilt.YLim,'color',[0,0.75,0]);
            end
            
            
            linkaxes([app.hAxRaw,app.hAxNorm,app.hAxFilt],'x')
            
            
            app.updateCurrentTrace();
        end
        
        function updateCurrentTrace(app)
            
            title(app.hAxRaw,['ROI ', num2str(app.tix)]);
            
            %raw traces
            app.hLraw(app.oldix).Color=app.lightgray;
            app.hLraw(app.oldix).LineWidth=app.bgLineWidth;
            app.hLraw(app.oldix).ZData=zeros(size(app.t));
            
            app.hLraw(app.tix).Color='k';
            app.hLraw(app.tix).LineWidth=app.boldLineWidth; %bold
            app.hLraw(app.tix).ZData=ones(size(app.t)); %on top of others
            
            %raw traces
            app.hLnorm(app.oldix).Color=app.lightgray;
            app.hLnorm(app.oldix).LineWidth=app.bgLineWidth;
            app.hLnorm(app.oldix).ZData=zeros(size(app.t));
            
            app.hLnorm(app.tix).Color='k';
            app.hLnorm(app.tix).LineWidth=app.boldLineWidth; %bold
            app.hLnorm(app.tix).ZData=ones(size(app.t)); %on top of others
            
            %filtered traces
            app.hLfilt(app.oldix).Color=app.lightgray;
            app.hLfilt(app.oldix).LineWidth=app.bgLineWidth;
            app.hLfilt(app.oldix).ZData=zeros(size(app.t));
            
            app.hLfilt(app.tix).Color='k';
            app.hLfilt(app.tix).LineWidth=app.boldLineWidth;
            app.hLfilt(app.tix).ZData=ones(size(app.t));
            
            delete(findobj(app.hFig,'Tag','spMarks'))
            
            hold(app.hAxRaw,'on')
            hold(app.hAxNorm,'on')
            hold(app.hAxFilt,'on')
            for interval=1:size(app.tIntervals,1)
                Pts=app.trace{interval}.points(app.tix);
                if ~isempty(Pts.tUp)
                    plot(app.hAxRaw,Pts.tUp,Pts.xUpR,'bs','Tag','spMarks');
                    plot(app.hAxNorm,Pts.tUp,Pts.xUpN,'bs','Tag','spMarks');
                    plot(app.hAxFilt,Pts.tUp,Pts.xUp,'bs','Tag','spMarks');
                    plot(app.hAxFilt,app.tIntervals(interval,:),Pts.thrUp*[1,1],'r--','Tag','spMarks');
                end
                if ~isempty(Pts.tDown)
                    plot(app.hAxRaw,Pts.tDown,Pts.xDownR,'bo','Tag','spMarks');
                    plot(app.hAxNorm,Pts.tDown,Pts.xDownN,'bo','Tag','spMarks');
                    plot(app.hAxFilt,Pts.tDown,Pts.xDown,'bo','Tag','spMarks');
                    if Pts.thrUp~=Pts.thrDown
                        plot(app.hAxFilt,app.tIntervals(interval,:),Pts.thrDown*[1,1],'b--','Tag','spMarks');
                    end
                end
                if ~isempty(Pts.tMax)
                    plot(app.hAxRaw,Pts.tMaxR,Pts.xMaxR,'rv','Tag','spMarks');
                    plot(app.hAxNorm,Pts.tMaxN,Pts.xMaxN,'rv','Tag','spMarks');
                    plot(app.hAxFilt,Pts.tMax,Pts.xMax,'rv','Tag','spMarks');
                    plot(app.hAxFilt,app.tIntervals(interval,:),Pts.thrUp*[1,1],'r--','Tag','spMarks');
                end
                if ~isempty(Pts.tMin)
                    plot(app.hAxRaw,Pts.tMinR,Pts.xMinR,'r^','Tag','spMarks');
                    plot(app.hAxNorm,Pts.tMinN,Pts.xMinN,'r^','Tag','spMarks');
                    plot(app.hAxFilt,Pts.tMin,Pts.xMin,'r^','Tag','spMarks');
                    plot(app.hAxFilt,app.tIntervals(interval,:),Pts.thrUp*[1,1],'r--','Tag','spMarks');
                end
                
            end
            hold(app.hAxRaw,'off')
            hold(app.hAxNorm,'off')
            hold(app.hAxFilt,'off')
            
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
                    disp('processing data...')
                    app.processData();
                    app.plotData();
                    
                case {'t'}
                    
                    
                case {'l'}
                    app.loadData();
                    app.processData();
                    app.plotData();
                    
                case {'s'}
                    disp('save placeholder')
                    
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
            
            if ~exist('filename','var')||isempty(filename)|| ~exist(filename,'file')
                [filename,path]=uigetfile({'*.xls*'});
                filename=[path,filename];
            end
            
            [app.datafilepath,app.datafilename]=fileparts(filename);
            
            [data,header]=xlsread(filename);
            time=data(:,1);
            X=data(:,2:end);
            
            %decode header info
            opts=[]; %placeholder
            
            %interpolate any missing data using neighboring time points
            rowIsNan=any(isnan(X),2);
            timeIsNan=isnan(time);
            
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
            
            
            if ~exist('tCond','var')
                tCond=[];
            end
            
            app.dt=mode(diff(app.t));
            app.fs=1/app.dt;
            %             tOmit=0;
            app.nOmit=round(app.tOmit/app.dt);
            tCond=[app.t(1),tCond(:)',app.t(end)];
            for i=1:length(tCond)-1
                app.tIntervals(i,:)=[tCond(i),tCond(i+1)];
            end
            
            app.nTraces=size(X,2);
            app.oldix=1;
            app.tix=1;
            
            app.hFig.Name=['Oscillation Analyzer - ',app.datafilename];
        end
        
        %save results/params: MAT + XLS? ask for save name
        function saveResults(app)
            %collect results into a table for writing to file
            
            
            %uiputfile
        end
    end
    
end