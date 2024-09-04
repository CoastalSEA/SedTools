classdef TideParams < muiPropertyUI                
%
%-------class help------------------------------------------------------===
% NAME
%   TideParams .m
% PURPOSE
%   Class for input parameters for the SimpleTide model
% USAGE
%   obj = TideParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI. Based on STparam from SimpleTide in ModelUI
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%      
    properties (Hidden)
        %abstract properties in muiPropertyUI to define input parameters
        PropertyLabels = {'Mean tide level (mOD)',...
                          'Tidal amplitude (m)',...
                          'Tidal period (hours)',...
                          'Elevation phase (deg)',...
                          'Velocity amplitude (m/s)',...
                          'Velocity phase (deg)',...
                          'M2 tidal amplitude (m)',...
                          'S2 tidal amplitude (m)',...
                          'O1 tidal amplitude (m)'};
        %abstract properties in muiPropertyUI for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        MSL0 = 0                  %mean tidel level to ordnance datum (mOD)
        TidalAmp                  %tidal elevation amplitude (m)
        TidalPeriod = 12.4        %tidal period (hours)
        ElevPhase = 0             %phase of elevation (ie k.x) (rads)
        VelocityAmp               %tidal velocity amplitude (m/s)
        VelocityPhase             %phase of velocity (ie k.x+phi) (rads)
        M2amplitude               %M2 tidal amplitude (m)
        S2amplitude               %S2 tidal amplitude (m)
        O1amplitude               %O1 tidal amplitude (m) 
    end    

    properties (Dependent)
        AngularFrequency          %tidal angular frequency (rads/s)
    end
%%   
    methods (Access=protected)
        function obj = TideParams(mobj)            
            %constructor code:            
            %values defined in UI function setTabProperties used to assign
            %the tabname and position on tab for the data to be displayed
            obj = setTabProps(obj,mobj);  %muiPropertyUI fcn
        end 
    end
%%  
    methods (Static)  
        function obj = setInput(mobj,editflag)
            %gui for user to set Parameter Input values
            classname = 'TideParams';           
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = TideParams(mobj);          
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

    methods
        function omega = get.AngularFrequency(obj)
            %tidal angular frequecy based on tidal period (rads/s)
            omega = 2*pi/(obj.TidalPeriod*3600);  
        end   
    end
end