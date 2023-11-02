classdef TransportAnalysis < muiDataSet                      
%
%-------class help---------------------------------------------------------
% NAME
%   TransportAnalysis.m
% PURPOSE
%   Class description - Class for Transport analysis in SedTools App
%
% SEE ALSO
%   muiDataSet
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%     
    properties
        %inherits Data, RunParam, MetaData and CaseIndex from muiDataSet
        %Additional properties:     
    end
    
    methods (Access = private)
        function obj = TransportAnalysis()                    
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
            obj = TransportAnalysis;                         
            dsp = modelDSproperties(obj);
            
            %now check that the input data has been entered
            if ~isValidModel(mobj, metaclass(obj).Name)  
                warndlg('Use Setup to define model input parameters');
                return;
            end
            muicat = mobj.Cases;
            %assign the run parameters to the model instance
            setRunParam(obj,mobj); 
            
%--------------------------------------------------------------------------
% Model code 
%--------------------------------------------------------------------------
            [results,meta,modeltime] = totaltransport_model(mobj);
%--------------------------------------------------------------------------
% Assign model output to a dstable using the defined dsproperties meta-data
%--------------------------------------------------------------------------                   
            %each variable should be an array in the 'results' cell array
            %if model returns single variable as array of doubles, use {results}
            dst = dstable(results{:},'RowNames',modeltime,'DSproperties',dsp);
%--------------------------------------------------------------------------
% Save results
%--------------------------------------------------------------------------                        
            %assign metadata about model
            dst.Source = metaclass(obj).Name;
            dst.MetaData = meta;
            %save results
            setDataSetRecord(obj,muicat,dst,'model');
            getdialog('Run complete');
        end
    end
%%
    methods
        function tabPlot(obj,src) %abstract class for muiDataSet
            %generate plot for display on Q-Plot tab
            
            %add code to define plot format or call default tabplot using:
            tabDefaultPlot(obj,src);
        end
    end 
%%    
    methods (Access = private)
        function dsp = modelDSproperties(~) 
            %define a dsproperties struct and add the model metadata
            dsp = struct('Variables',[],'Row',[],'Dimensions',[]); 
            %define each variable to be included in the data table and any
            %information about the dimensions. dstable Row and Dimensions can
            %accept most data types but the values in each vector must be unique
            
            %struct entries are cell arrays and can be column or row vectors
            dsp.Variables = struct(...                       % <<Edit metadata to suit model
                'Name',{'Depth','Velocity','Taux','Taur','Qt'},...
                'Description',{'Depth','Velocity','Maximum bed shear stress',...
                               'RMS bed shear stress','Total transport'},...
                'Unit',{'m','m/s','N/m^2','N/m^2','m^2/s'},...
                'Label',{'Depth (m)','Velocity (m/s)','Bed shear stress (N/m^2)',...
                         'Bed shear stress (N/m^2)','Total transport (m^3/m-width/s)'},...
                'QCflag',{'model','model','model','model','model'}); 
            dsp.Row = struct(...
                'Name',{'Time'},...
                'Description',{'Time'},...
                'Unit',{'h'},...
                'Label',{'Time (d)'},...
                'Format',{'d'});        
            dsp.Dimensions = struct(...    
                'Name',{''},...
                'Description',{''},...
                'Unit',{''},...
                'Label',{''},...
                'Format',{''});  
        end
    end    
end