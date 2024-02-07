classdef SedTools < muiModelUI
%
%-------class help---------------------------------------------------------
% NAME
%   SedTools.m
% PURPOSE
%   Main UI for SedTools implementation using ModelUI.
% NOTES
%   based on cal_settling by Urs Neumeier, 2005, 
%   http://neumeier.perso.ch/matlab/cal_settling.html
% SEE ALSO
%   Abstract class muiModelUI.m, ModelUI.m and tools provided in muitoolbox
%
% Author: Ian Townend
% CoastalSEA (c) Jan 2021
%--------------------------------------------------------------------------
% 
    properties  (Access = protected)
        %implement properties defined as Abstract in muiModelUI
        vNumber = '2.0'
        vDate   = 'Nov 2023'
        modelName = 'SedTools'  
        %Properties defined in muiModelUI that need to be defined in setGui
        % ModelInputs  %classes required by model: used in isValidModel check 
        % DataUItabs   %struct to define type of muiDataUI tabs for each use                         
    end
    
    methods (Static)
        function obj = SedTools
            %constructor function initialises GUI  
            isok = check4muitoolbox(obj);
            if ~isok, return; end
            %
            obj = setMUI(obj);  
        end
    end
%% ------------------------------------------------------------------------
% Definition of GUI Settings
%--------------------------------------------------------------------------  
    methods (Access = protected)
        function obj = setMUI(obj)
            %initialise standard figure and menus   
            %classes required to run model; format:          
            %obj.ModelInputs.<model classname> = {'Param_class1',Param_class2',etc}
            obj.ModelInputs.SettlingAnalysis = {'SettlingParams'};    
            obj.ModelInputs.TransportAnalysis = {'TransportParams','TideParams','SiteParams','RunParams'};    
            %tabs to include in DataUIs for plotting and statistical analysis
            %select which of the options are needed and delete the rest
            %Plot options: '2D','3D','4D','2DT','3DT','4DT'
            obj.DataUItabs.Plot = {'2D','3D'};  
            %Statistics options: 'General','Timeseries','Taylor','Intervals'
            obj.DataUItabs.Stats = {'General','Taylor'}; 
            modelLogo = 'SedTools_logo.jpg'; %default splash figure 
            initialiseUI(obj,modelLogo); %initialise menus and tabs   
        end    
%% ------------------------------------------------------------------------
% Definition of Menu Settings
%--------------------------------------------------------------------------
        function menu = setMenus(obj)
            %define top level menu items and any submenus
            %MenuLabels can any text but should avoid these case-sensitive 
            %reserved words: "default", "remove", and "factory". If label 
            %is not a valid Matlab field name this the struct entry
            %is modified to a valid name (eg removes space if two words).
            %The 'gcbo:' Callback text triggers an additional level in the 
            %menu. Main menu labels are defined in sequential order and 
            %submenus in order following each brach to the lowest level 
            %before defining the next branch.         
                                                              % << Edit menu to suit model 
            MenuLabels = {'File','Tools','Project','Setup','Run',...
                                                        'Analysis','Help'};
            menu = menuStruct(obj,MenuLabels);  %create empty menu struct
            %
            %% File menu --------------------------------------------------
             %list as per muiModelUI.fileMenuOptions
            menu.File.List = {'New','Open','Save','Save as','Exit'};
            menu.File.Callback = repmat({@obj.fileMenuOptions},[1,5]);
            
            %% Tools menu -------------------------------------------------
            %list as per muiModelUI.toolsMenuOptions
            menu.Tools(1).List = {'Refresh','Clear all'};
            menu.Tools(1).Callback = {@obj.refresh, 'gcbo;'};  
            
            % submenu for 'Clear all'
            menu.Tools(2).List = {'Model','Figures','Cases'};
            menu.Tools(2).Callback = repmat({@obj.toolsMenuOptions},[1,3]);

            %% Project menu -----------------------------------------------
            menu.Project(1).List = {'Project Info','Cases','Export/Import'};
            menu.Project(1).Callback = {@obj.editProjectInfo,'gcbo;','gcbo;'};
            
            %list as per muiModelUI.projectMenuOptions
            % submenu for Scenarios
            menu.Project(2).List = {'Edit Description','Edit Data Set',...
                                    'Save Data Set','Delete Case','Reload Case',...
                                    'View Case Settings'};                                               
            menu.Project(2).Callback = repmat({@obj.projectMenuOptions},[1,6]);
            
            % submenu for 'Export/Import'                                          
            menu.Project(3).List = {'Export Case','Import Case'};
            menu.Project(3).Callback = repmat({@obj.projectMenuOptions},[1,2]);
            
            %% Setup menu -------------------------------------------------
            menu.Setup(1).List = {'Settling Parameters','Transport Parameters','Model Constants'};                                                                      
            menu.Setup(1).Callback = {@obj.setupMenuOptions,'gcbo;',@obj.setupMenuOptions};
            menu.Setup(1).Separator = {'off','off','on'};%separator preceeds item

            menu.Setup(2).List = {'Tidal Parameters','Sediment Parameters',...
                                        'Site Parameters','Run Parameters'};
            menu.Setup(2).Callback = repmat({@obj.setupMenuOptions},[1,4]);
            %% Run menu ---------------------------------------------------
            menu.Run(1).List = {'Settling Analysis','Transport Analysis',...
                                    'Sediment Properties','Derive Output'};
            menu.Run(1).Callback = repmat({@obj.runMenuOptions},[1,4]);
            menu.Run(1).Separator = {'off','off','off','on'};%separator preceeds item

            %% Plot menu --------------------------------------------------  
            menu.Analysis(1).List = {'Plots','Statistics'};
            menu.Analysis(1).Callback = repmat({@obj.analysisMenuOptions},[1,2]);
            
            %% Help menu --------------------------------------------------
            menu.Help.List = {'Documentation','Manual'};
            menu.Help.Callback = repmat({@obj.Help},[1,2]);
            
        end
        
%% ------------------------------------------------------------------------
% Definition of Tab Settings
%--------------------------------------------------------------------------
        function [tabs,subtabs] = setTabs(obj)
            %define main tabs and any subtabs required. struct field is 
            %used to set the uitab Tag (prefixed with sub for subtabs). 
            %Order of assignment to struct determines order of tabs in figure.
            %format for tabs: 
            %    tabs.<tagname> = {<tab label>,<callback function>};
            %format for subtabs: 
            %    subtabs.<tagname>(i,:) = {<subtab label>,<callback function>};
            %where <tagname> is the struct fieldname for the top level tab.
            tabs.Cases  = {'   Cases  ',@obj.refresh};      
            tabs.Inputs = {'  Inputs  ',''};
            subtabs.Inputs(1,:) = {' Settling ',@obj.InputTabSummary};
            subtabs.Inputs(2,:) = {' Transport ',@obj.InputTabSummary};
            tabs.Plot   = {'  Q-Plot  ',@obj.getTabData};
            tabs.Stats = {'   Stats   ',@obj.setTabAction};
        end 
 %%
        function props = setTabProperties(~)
            %define the tab and position to display class data tables
            %props format: {class name, tab tag name, position, ...
            %               column width, table title}                                                    
            props = {...                                     
               'SettlingParams','Settling',[0.95,0.5],{180,60},'Settling parameters:';...
               'TideParams','Transport',[0.95,0.45],{180,60},'Tidal parameters:';...
               'TransportParams','Transport',[0.95,0.95],{180,60},'Transport parameters:';...               
               'SiteParams','Transport',[0.5,0.95],{180,60},'Site parameters:';...
               'RunParams','Transport',[0.4,0.45],{180,60},'Run parameters:'};
        end      
 %%
        function setTabAction(obj,src,cobj)
            %function required by muiModelUI and sets action for selected
            %tab (src)
            msg = 'No results to display';
            switch src.Tag                   
                case 'Plot' 
                    tabPlot(cobj,src);
                case 'Stats'
                    lobj = getClassObj(obj,'mUI','Stats',msg);
                    if isempty(lobj), return; end
                    tabStats(lobj,src);    
            end
        end
        
%% ------------------------------------------------------------------------
% Callback functions used by menus and tabs
%-------------------------------------------------------------------------- 
        %% File menu ------------------------------------------------------
        %use default menu functions defined in muiModelUI
            
        %% Tools menu -----------------------------------------------------
        %use default menu functions defined in muiModelUI
                
        %% Project menu ---------------------------------------------------
        %use default menu functions defined in muiModelUI   
        %% Setup menu -----------------------------------------------------
        function setupMenuOptions(obj,src,~)
            %callback functions for data input
            tabname = 'Transport';
            switch src.Text
                case 'Settling Parameters'                
                    SettlingParams.setInput(obj);  
                    tabname = 'Settling';
                case 'Tidal Parameters'
                    TideParams.setInput(obj);                   
                case 'Sediment Parameters'
                    TransportParams.setInput(obj); 
                case 'Site Parameters'  
                    SiteParams.setInput(obj); 
                case 'Run Parameters'  
                    RunParams.setInput(obj);     
                case 'Model Constants'
                    obj.Constants = setInput(obj.Constants);
                    tabname = [];
            end
            %
            tabUpdate(obj,tabname);
        end             
        %% Run menu -------------------------------------------------------
        function runMenuOptions(obj,src,~)
            %callback functions to run model
            switch src.Text                   
                case 'Settling Analysis'                       
                    SettlingAnalysis.runModel(obj); 
                case 'Transport Analysis'
                    TransportAnalysis.runModel(obj);
                case 'Sediment Properties'
                    sediment_properties();
                case 'Derive Output'
                    obj.mUI.Manip = muiManipUI.getManipUI(obj);
            end 
            DrawMap(obj);
        end 
        %% Analysis menu ------------------------------------------------------
        function analysisMenuOptions(obj,src,~)
            switch src.Text
                case 'Plots'
                    obj.mUI.PlotsUI = muiPlotsUI.getPlotsUI(obj);
                case 'Statistics'
                    obj.mUI.StatsUI = muiStatsUI.getStatsUI(obj);
            end            
        end

        %% Help menu ------------------------------------------------------    
        function Help(~,src,~)
            %menu to access online documentation and manual pdf file
            switch src.Text
                case 'Documentation'
                    doc sedtools   %must be name of html help file  
                case 'Manual'
                    sdt_open_manual;
            end
        end 
        %% Check that toolboxes are installed------------------------------
        function isok = check4muitoolbox(~)
            %check that dstoolbox and muitoolbox have been installed
            fname = 'dstable.m';
            dstbx = which(fname);
        
            fname = 'muiModelUI.m';
            muitbx = which(fname);
        
            if isempty(dstbx) && ~isempty(muitbx)
                warndlg('dstoolbox has not been installed')
                isok = false;
            elseif ~isempty(dstbx) && isempty(muitbx)
                warndlg('muitoolbox has not been installed')
                isok = false;
            elseif isempty(dstbx) && isempty(muitbx)
                warndlg('dstoolbox and muitoolbox have not been installed')
                isok = false;
            else
                isok = true;
            end
        end        
    end
%%
%%
    methods (Access=private)
        function tabUpdate(obj,tabname)
            %update tab used for properties if required
            if ~isempty(tabname)
                tabsrc = findobj(obj.mUI.Tabs,'Tag',tabname);
                InputTabSummary(obj,tabsrc);
            end 
        end
    end
end