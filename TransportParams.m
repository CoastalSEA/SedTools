classdef TransportParams < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   TransportParams.m
% PURPOSE
%   Class for total transport sediment parameters
% USAGE
%   obj = TransportParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%    
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Temperature (degC)',...
                          'Salinity (ppt)',...
                          'Median grain size, D50 (m)',...
                          'D90 grain size (m)',...
                          'Percentage mud content (%)',...
                          'Bulk density of seabed (kg/m3)'}
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        WaterTemp = 10    %Temperature (degC)
        Salinity = 32     %Salinity (ppt)'
        D50               %Median grain size, D50 (m)
        D90               %D90 grain size (m)
        PercentMud = 0    %Percentage mud content (%)
        BedDensity        %Bulk density of seabed (kg/m3)                                            
    end      

%%   
    methods (Access=protected)
        function obj = TransportParams(mobj)  %instantiate class
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
            classname = 'TransportParams';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = TransportParams(mobj);             
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