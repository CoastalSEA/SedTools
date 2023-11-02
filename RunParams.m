classdef RunParams < muiPropertyUI              
%
%-------class help---------------------------------------------------------
% NAME
%   RunParams.m
% PURPOSE
%   Class for to set run time properties - used in SedTools
% USAGE
%   obj = RunParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {' Duration of simulation (days)',...
                          ' Time step (mins)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        Duration = 14.77           %duration (days)
        Tinterval = 10             %time interval for simulation (mins)
    end    

%%   
    methods (Access=protected)
        function obj = RunParams(mobj)         
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI function
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'RunParams';           
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = RunParams(mobj);             
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
end