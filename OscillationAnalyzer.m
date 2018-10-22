classdef OscillationAnalyzer<handle
    %OscillationAnalyzer is a graphical interface for performing
    %processing and feature extraction on timeseries data.
    
    
    %data and processing properties
    properties
        
        expt
        
        featmethod='threshold'
        thrFrac=[0.55,0.45]
        delta=0.25
        Tbig=15
        Tsmall=4
        
        xfeat
        yfeat
        
    end
    
    %GUI properties
    properties
        
        hFig
        axRaw
        axNorm
        axDT
        axFilt
        axPSD
        axFeat
        
        tix=1 %focus trace
        oldix=1
        lightgray=[0.75,0.75,0.75]
        bgLineWidth=0.5
        boldLineWidth=1.5
        
        %uicontrols
        
        
    end
    
    %constructor
    methods
        function app=OscillationAnalyzer(filename,condNames,condTimes)
            %assume data is either filename, or data matrix [t,X], columns
            %are observations.
            
            app.expt=Experiment(filename);
            
            app.expt.defineSegments(condNames,condTimes);
            
            app.buildGUI(); %set up display elements
            
            app.processData();
            
            app.xfeat='Tpsd';
            app.yfeat='Tmean';
%             app.xfeat='segment';
%             app.yfeat=app.expt.fnames{1};
            
            app.plotData(); %populate the axes
        end
    end
    
    %data processing and feature detection functions
    methods
        
        
        function processData(app)
            
            app.expt.normalize('devmean2',[],0)

            persegment=false;
            app.expt.detrend('sgolay',app.Tbig,persegment,0)

            app.expt.filter('gaussian',app.Tsmall,0)
            
            switch app.featmethod
                case 'peaks'
                    app.expt.compute_features('peaks',app.delta);
                case 'threshold'
                    app.expt.compute_features('threshold',app.thrFrac);
            end
            
            
        end
        
    end
    
    %callbacks and other GUI functions
    methods
        %modify default analysis parameters
        %TODO: UI table? pops up if needed?
        function setParams(app)
        end
        
        
        %set up all the graphical components of the UI
        function buildGUI(app)
            app.hFig=figure('Name',['Oscillation Analyzer - ',app.expt.name],'NumberTitle','off');
            app.hFig.KeyPressFcn=@app.keyPressDecoder;
            pos=get(groot,'screensize');
            pos=[50,50,pos(3:4)-150];
            app.hFig.Position=pos;
            app.hFig.Units='normalized';
            
            xl=0.075; yl=0.075; wl=0.5; hl=0.2; gl=0.025;
            xr=0.65; wr=0.3; hr=0.3; gr=0.05;
            app.axRaw=axes('Position',[xl,yl+3*(hl+gl),wl,hl]);
            app.axNorm=axes('Position',[xl,yl+2*(hl+gl),wl,hl]);
            app.axDT=axes('Position',[xl,yl+hl+gl,wl,hl]);
            app.axFilt=axes('Position',[xl,yl,wl,hl]);
            app.axPSD=axes('Position',[xr,yl+hr+gr,wr,hr]);
            app.axFeat=axes('Position',[xr,yl,wr,hr]);
            
            
            
        end
        
        function plotData(app)
            showPts=true;
            app.expt.plotTrace(app.axRaw,'raw',showPts,app.tix)
            app.expt.plotTrace(app.axNorm,'norm',showPts,app.tix)
            app.expt.plotTrace(app.axDT,'detrend',showPts,app.tix)
            app.expt.plotTrace(app.axFilt,'filt',showPts,app.tix)
            
            set([app.axRaw,app.axNorm,app.axDT],'Xticklabel','');
            
            linkaxes([app.axRaw,app.axNorm,app.axDT,app.axFilt],'x')
            
            app.expt.plotPeriodogram(app.axPSD,app.tix)
            
            app.expt.plotFeatures(app.axFeat,app.xfeat,app.yfeat,app.tix)
            
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
                        app.plotData();
                    end
                case {'rightarrow'}
                    if app.tix<app.expt.nX
                        app.oldix=app.tix;
                        app.tix=app.tix+1;
                        app.plotData();
                    end
                case {'p'}
                    disp('processing data...')
                    app.processData();
                    app.plotData();
                    
                case {'t'}
                    prompt = {'Enter trend window duration:','Enter smoothing window duration:'};
                    title = 'Input parameters';
                    dims = [1 35];
                    definput = {num2str(app.Tbig),num2str(app.Tsmall)};
                    answer = inputdlg(prompt,title,dims,definput);
                    app.Tbig=str2double(answer{1});
                    app.Tsmall=str2double(answer{2});
                    app.processData();
                    app.plotData();
                    
                case {'l'}
                    app.loadData();
                    app.processData();
                    app.plotData();
                    
                case {'s'}
                    disp('save placeholder')
                    
                case '1'
                    app.featmethod='peaks';
                    
                case '2'
                    app.featmethod='threshold';
                    
                    
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
            
            app.expt=Experiment(filename);
            
            app.hFig.Name=['Oscillation Analyzer - ',app.datafilename];
        end
        
        %save results/params: MAT + XLS? ask for save name
        function saveResults(app)
            %collect results into a table for writing to file
            
            
            %uiputfile
        end
    end
    
end