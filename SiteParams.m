classdef SiteParams < muiPropertyUI
%
%-------class help---------------------------------------------------------
% NAME
%   SiteParams.m
% PURPOSE
%   Class for total transport site parameters
% USAGE
%   obj = SiteParams.setInput(mobj); %mobj is a handle to Main UI
% SEE ALSO
%   inherits muiPropertyUI
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%    
    properties (Hidden)
        %abstract properties in PropertyInterface to define input variables
        PropertyLabels = {'Significant wave height (m)',...
                          'Peak wave period (s)',...
                          'Angle between flow and waves (deg)',...
                          'Bed elevation (mAD)',...
                          'Bed slope (1:bs. +ve if flow is upslope)',...
                          'Bed roughness coefficient (-)',...                          
                          'Ripple height (m)',...
                          'Ripple length (m)'}
        %abstract properties in PropertyInterface for tab display
        TabDisplay   %structure defines how the property table is displayed 
    end
    
    properties
        Hs                    %Significant wave height (m)
        Tp                    %Peak wave period (s)
        Phid                  %Angle between flow and waves (deg)
        zBed                  %Bed elevation (mAD) 
        BedSlope              %Bed slope(1:bs) in streamwise direction (+ve if flow is upslope)
        BedRoughness = 0      %Bed roughness coefficient (-)          
        RippleHeight = 0      %Ripple height - adjusts d50 and d90 if >0 (m)
        RippleLength = 0      %Ripple length (m)
    end   

    properties (Dependent)
        Depth2mtl             %water depth relavtive to mean tide level
    end

%%   
    methods (Access=protected)
        function obj = SiteParams(mobj)  %instantiate class
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
            classname = 'SiteParams';               
            obj = getClassObj(mobj,'Inputs',classname);
            if isempty(obj)
                obj = SiteParams(mobj);             
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
        function depth = get.Depth2mtl(obj)
            %use elevation of mtl and bed elevation to get depth
            depth = obj.MSL0-obj.zBed;
        end
    end
end