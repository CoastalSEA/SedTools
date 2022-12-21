%% SedTools
% Matlab App for analysis of settling column data

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

%% SedTools classes
% * *SedTools* - defines the behaviour of the main UI.
% * *SettlingParams* - defines the input parameters.
% * *SettlingAnalysis* - handles the running of settling_column function,
% saving of results and display on the Q-Plot tab.

%% SedTools functions
% *settling_column* - convert the record from the settling columns to 
% grain size in units of mm or phi.

%% Manual
% The <matlab:sdt_open_manual manual> provides further details of setup and 
% configuration of the model. Sample input files can be found in
% the example folder <matlab:sdt_example_folder here>. 

%% See Also
% <matlab:doc('muitoolbox') muitoolbox>, <matlab:doc('dstoolbox') dstoolbox>.
	