%% SedTools
% Matlab App for analysis of settling column data and undertaking various sediment property and transport analyses.

%% Licence
% The code is provided as Open Source code (issued under a GNU General 
% Public License).

%% Requirements
% SedTools is written in Matlab(TM) and requires v2016b, or later. In addition, 
% SedTools requires both the <matlab:doc('dstoolbox') dstoolbox> and the 
% <matlab:doc('muitoolbox') muitoolbox>.

%% Background
% Utility to analyse settling column data based on the analysis method
% developed by Urs Neumeier, 2005. The source code of the core program and 
% further details see <http://neumeier.perso.ch/matlab/cal_settling.html web site>. 
% The total load transport model is based on Sousby & van_Rijn in Dynamics
% of Marine Sands, 1997.

%% SedTools classes
% * *SedTools* - defines the behaviour of the main UI.
% * *SettlingParams* - defines the input parameters.
% * *SettlingAnalysis* - handles the running of settling_column function,
% saving of results and display on the Q-Plot tab.
% * *RunParams* – run time parameters for total load transport model.
% * *SiteParams* – wave and bed conditions at the site being modelled.
% * *TideParams* – data entry of tidal properties for use in the tidal emulator model.
% * *TransportAnalysis* – class to handle and call the total transport model.
% * *TransportParams* – sediment parameters used in total load transport model

%% SedTools functions
% * *settling_column* - convert the record from the settling columns to 
% grain size in units of mm or phi.
% * *sandwave* – compute the height and wave length of bed ripples using van Rijn (1984)
% * *svr_transport* – function to calculate the total load transport using the equations of Soulsby-van Rijn in Dynamics of Marine Sands.
% * *totaltransport_model* – run the total load transport model for SedTools App
% * *transport_plot* –  plot results from the total load transport model in the SedTools App

%% muiAppCoastalFcns
% * *celerity* – calculate the wave celerity using Hunt's equation
% * *fluidprops* – calculate fluid density and kinematic viscosity based on salinity and temperature.
% * *hb_break* – wave height after breaking for given water depth d
% * *hs_break* – significant wave height after breaking for given water depth d
% * *sedprops* –returns one of a range of sediment properties based on user selection 
% * *sediment_properties* – calculate a set of sediment properties based on bulk density
% * *simple_tide* – compute a tidal water level time series using main constituents scaled to required tidal amplitude  
% * *tau_bed* – compute bed shear stress under combined wave-current action
% * *tau_crit* – calculate the critical erosion shear stress and erosion rate for sand, mud or mixed sediments

%% Manual
% The <matlab:sdt_open_manual manual> provides further details of setup and 
% configuration of the model. Sample input files can be found in
% the example folder <matlab:sdt_example_folder here>. 

%% See Also
% <matlab:doc('muitoolbox') muitoolbox>, <matlab:doc('dstoolbox') dstoolbox>.
	