classdef SettlingAnalysis < muiDataSet       
%
%-------class help---------------------------------------------------------
% NAME
%   SettlingAnalysis.m
% PURPOSE
%   Class description - Class for Settling Analysis run as a muitoolbox app
%
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:   
        Stats    %table for statistical summary of results
    end
    
    methods (Access = private)
        function obj = SettlingAnalysis()                      
            %class constructor
        end
    end      
%%
    methods (Static)        
%--------------------------------------------------------------------------
% Model implementation
%--------------------------------------------------------------------------         
        function obj = runModel(mobj)
            %function to run a simple 2D diffusion model
            obj = SettlingAnalysis;                        
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            %isValidModel checks the InputHandles defined in ModelUI
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            %may need to be after input data selection to capture caserecs
            setRunParam(obj,mobj); 
            %--------------------------------------------------------------
            % Model code 
            %--------------------------------------------------------------
            inp = mobj.Inputs.SettlingParams;   
            dsp.Dimensions.Unit = inp.OutputUnits;
            %get files of data to be analysed
            [fname,path,nfiles] = getfiles('MultiSelect','on',...
                                     'FileType',{'*.dat;*.txt'},...
                                    'PromptText','Select sediment files' );
            %to average the results from several measurements (files) 
            %data is loaded as multiple columns.
            data{1,nfiles} = [];
            for i=1:nfiles
                filename = [path,fname{i}];
                data(i) = getData(obj,filename);    %get data returns cell            
            end
            [grainsize,results,stats,metatxt]=settling_column(data,inp);
            %--------------------------------------------------------------
            % Assign output to a dstable using defined dsproperties meta-data
            %--------------------------------------------------------------               
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            results = cellfun(@transpose,results,'UniformOutput',false);
            dst = dstable(results{:},'DSproperties',dsp);
            dst.Dimensions.D = grainsize{1};     %grain size
            %--------------------------------------------------------------
            % Save results
            %--------------------------------------------------------------                      
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = metatxt;
            obj.Stats = cell2table(stats(:,2),'RowNames',stats(:,1),...
                'VariableNames',{'Value'});
            %save results
            casedesc = regexp(fname{1},'[.]','split');
            setDataSetRecord(obj,muicat,dst,'analysis',casedesc(1));
            getdialog('Run complete');
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            if strcmp(src.Tag,'FigButton')
                hfig = figure('Tag','PlotFig');
                ax = axes('Parent',hfig,'Tag','PlotFig','Units','normalized');
                grainSizePlot(obj,ax);
            else
                ht = findobj(src,'Type','panel');
                delete(ht);
                hp1 = uipanel('Parent',src,'Position',[0,0,0.7,1],'Tag','PlotTab');
                ax = axes('Parent',hp1,'Tag','MyPlot');
                grainSizePlot(obj,ax);
                
                %now display table of stats
                stats = obj.Stats; 
                hp2 = uipanel('Parent',src,'Position',[0.7,0,0.3,1],...
                    'BorderType','none','Tag','PlotTab');
                ht = uitable(hp2,'RowName',stats.Properties.RowNames,...
                    'ColumnName',{'Value'},...
                    'ColumnWidth',{60},'Data',stats.Value,...
                    'Units','normalized','Position',[0.02,0,0.98,0.9]);
                ht.Position(3) = ht.Extent(3)+0.1;
                ht.Position(2) = 0.1;
                %button to copy data to clipboard
                uicontrol('Parent',hp2,'Style','pushbutton',...
                    'String','Copy to clipboard','UserData',stats,...
                    'TooltipString','Copy stats table to clipboard',...
                    'Units','normalized','Position',[0.2 0.02 0.6 0.06],...
                    'Callback',@copydata2clip);
                %button to create plot as stand-alone figure
                uicontrol('Parent',hp1,'Style','pushbutton',...
                    'String','>Figure','Tag','FigButton',...
                    'TooltipString','Create plot as stand alone figure',...
                    'Units','normalized','Position',[0.88 0.95 0.10 0.044],...
                    'Callback',@(src,evdat)tabPlot(obj,src));
            end
        end
    end 
%%    
    methods (Access = private)
        function data = getData(~,filename)
            %read wave data (read format is file specific).
            fid = fopen(filename, 'r');
            if fid<0
                errordlg('Could not open file for reading','File write error','modal')
                data = [];
                return;
            end
            %read numeric data            
            data = textscan(fid,'%f');
            if isempty(data)
                warndlg('No data. Check file format selected')
            end
            fclose(fid);
        end
%%
        function sizeclasses = getSizeClasses(~,units,grainsize)
            %recover class boundaries from size distribution central values
            if strcmp(units,'mm')
                grainsize = -log2(grainsize);                
            end
            gsint = grainsize(1)-grainsize(2);
            sizeclasses(:,1) = grainsize+gsint/2;
            sizeclasses(end+1,1) = grainsize(end)-gsint/2;
            if strcmp(units,'mm')                
                sizeclasses = 2.^(-sizeclasses); 
            end
        end        
%%
        function grainSizePlot(obj,ax)
            %default tab plot of grain-size distribution
            dst = obj.Data.Dataset;        
            %--------------------------------------------------------------
            %extract data vectors, eg:
            gs = dst.Dimensions.D;
            cf = dst.cumfreq;
            sf = dst.sedfreq;
            units = dst.DimensionUnits{1};
            szclass = getSizeClasses(obj,units,gs);
            
            %surface plot of 3D form
            if strcmp(units,'phi')
                ax.XDir = 'reverse';
                szclass = flipud(szclass);
                sf = flipud(sf);
            else
                ax.XScale = 'log';
            end
            xlabel(sprintf('%s (%s)',dst.DimensionLabels{1},units));
            yyaxis left
            %bar(ax,gs,sf);
            histogram(ax,'BinEdges',szclass','BinCounts',sf)
            ylabel('Grain-size frequency (%)');
            yyaxis right
            plot(ax,gs,cf,'-r','LineWidth',1);
            ylabel('Cummulative frequency distribution (%)');
            ylim([0,100]);
            title (dst.Description);           
            ax.Color = [0.96,0.96,0.96];  %needs to be set after plot 
        end
%%
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                       % <<Edit metadata to suit model
                'Name',{'cumfreq','sedfreq'},...
                'Description',{'Cummulative frequency','Grain-size frequency'},...
                'Unit',{'%','%'},...
                'Label',{'Cummulative frequency distribution','Grain-size frequency (%)'},...
                'QCflag',{'analysis','analysis'}); 
            dsp.Row = struct(...
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});        
            dsp.Dimensions = struct(...    
                'Name',{'D'},...
                'Description',{'Grain-size'},...
                'Unit',{''},...
                'Label',{'Grain-size'},...
                'Format',{''});  
        end
    end    
end