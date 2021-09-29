classdef SettlingParams < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   PropsInput_template.m
% PURPOSE
%   Class for settling model input parameters
% USAGE
%   obj = SettlingParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%    
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Temperature (°C)',...
                          'Sediment density (kg/m3)',...
                          'Column height (m)',...
                          'Recorder time-step (s)',...
                          'Period from start to ignore (s)',...
                          'Settling velocity method (1-4)',...
                          'Resolution for results output (1,2,4,8)',...
                          'Output units in ''phi'' or ''mm'''}
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        Temperature                %Temperature (°C)
        SedimentDensity = 2650     %Sediment density (default 2650 kg/m3)
        ColumnHeight = 1.768       %Height (metres, default 1.768 m)
        RecTimeStep = 0.1591       %Time-step of the recorder in seconds
        OmitPeriod = 0             %period to ignore from the start of the record (s)
        SettlingMethod = 4         %method used to compute settling velocity
                                   %1=Gibbs et al. (1971); 2=Baba and Komar (1981)
                                   %3=van Rijn (1993); 4=Soulsby (1997)
        Resolution = 8             %resolution for the results (fractions of phi)
        OutputUnits = 'phi'        %output units: phi or mm                                 
    end    

%%   
    methods (Access=protected)
        function obj = SettlingParams(mobj)  %instantiate class
            %constructor code: model specific properties for class
            %uses the values defined in UI function setTabProperties
            %inputs: obj-class handle; mobj-UI handle; <class name>
            obj = setTabProps(obj,mobj);        %PropertyInterface fcn
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'SettlingParams';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = SettlingParams(mobj);             
            end
            %use muiPropertyUI function to generate UI
            if nargin<2 || editflag
                %add nrec to limit length of props UI (default=12)
                obj = editProperties(obj);  
                %add any additional manipulation of the input here
            end
            setClassObj(mobj,'Inputs',classname,obj);
        end        
    end
%%        
        %add other functions to operate on properties as required
        
end