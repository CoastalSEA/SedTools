%% Menu Options
% Summary of the options available for each drop down menu.

%% File
% * *New*: clears any existing model (prompting to save if not already saved) and a popup dialog box prompts for Project name and Date (default is current date). 
% * *Open*: existing Asmita models are saved as *.mat files. User selects a model from dialog box.
% * *Save*: save a file that has already been saved.
% * *Save as*: save a file with a new or different name.
% * *Exit*: exit the program. The close window button has the same effect.

%% Tools
% * *Refresh*: updates Cases tab.
% * *Clear all > Project*: deletes the current project, including all Setup data and all Cases.
% * *Clear all > Figures*: deletes all results plot figures (useful if a large number of plots have been produced).
% * *Clear all > Cases*: deletes all Cases listed on the Cases tab but does not affect the model setup.

%% Project
% * *Project Info*: edit the Project name and Date
% * *Cases > Edit Description*: user selects a Case to edit the Case description.
% * *Cases > Edit DS properties*: initialises the  UI for editing Data Set properties (DSproperties).
% * *Cases > Edit Data Set*: initialises the Edit Data UI for editing data sets.
% * *Cases > Save*: user selects a data set to be saved from a list box of Cases and the is then prompted to name the file. The data are written to an Excel spreadsheet. 
% * *Cases > Delete*: user selects Case(s) to be deleted from a list box of Cases and results are then deleted (model setup is not changed).
% * *Cases > Reload*: user selects a Case to reload as the current parameter settings.
% * *Cases > View settings*: user selects a Case to display a table listing the parameters used for the selected Case. 
% * *Export/Import > Export*: user selects a Case class instance to export as a mat file.
% * *Export/Import > Import*: user selects an exported Case class instance (mat file) to be loaded.
%%
% *NB*: to export the data from a Case for use in another application 
% (eg text file, Excel, etc), use the *Project>Cases>Edit Data Set* option 
% to make a selection and then use the ‘Copy to Clipboard’ button to paste 
% the selection to the clipboard.

%% Setup
% * *Settling Parameters*: dialogue to define settling analysis input 
% parameters.
% * *Model Constants*: a number of constants are used in the model. Generally, the default values are appropriate but these can be adjusted and saved with the project if required.
%%
% * *Settling Parameters* - the inputs required are as follows:
%%
%     Temperature during the settling experiment (*).
%     Density of the sediment used (*).
%     Height of the settling column (*).
%     The time step used to record the settling data
%     The period from the start of the recording to omit from the analysis
%     Computation method for settling velocity (+)
%     The resolution to be used in the output (values are 1/phi)
%     Units to display output
%%
%     (*) If multiple data sets are being averaged, multiple values can be 
%     entered, with the number of values equal to the number of input files. 
%     Otherwise the first value is replicated for all data sets. 
%%
%     (+) Methods: 1 = Gibbs et al. (1971); 2 = Baba and Komar (1981); 3 = van Rijn (1993); 4 = Soulsby (1997). 
%%
%     *References*
%%
%     1. Gibbs R.J., Matthews M.D. and Link D.A. (1971) The relationship between sphere size and settling velocity. Journal of sedimentary Petrology, 41/1, 7-18.
%     2. Baba J. & Komar P.D. (1981) Measurements and analysis of settling velocities of natural quartz sand grains. Journal of sedimentary Petrolology, 51, 631-640.
%     3. van Rijn L.C. (1993) Principles of sediment transport in rivers, estuaries and coastal seas. Aqua Publications.
%     4. Soulsby R. (1997) Dynamics of marine sands: a manual for practical applications. Thomas Telford.
%%
% * *Transport Parameters* using the following sub menus:
%%
% * _Tidal Parameters_: define the tidal conditions, i.e., tidal amplitude,
% velocity, phase, etc.
% * _Sediment Parameters_: define sediment properties used in the model,
% i.e. D50, D90, density of bed, etc. Note that grain sizes are in metres.
% * _Site Parameters_: define the conditions at the site being examined,
% i.e., wave conditions, bed description, etc. 
% * _Run Parameters_: specify the duration and time to be used, i.e.,
% duration and time step.

%% Run
% * *Settling Analysis*: runs analysis of settlingm column output files.
% Multiple files can be selected to take the analyse the average of across
% all selected files. On completion the first file name is given as the
% default Case description, and can be edited before being saved. The 
% results are added to the listing on the Cases tab and can viewed on the Q-Plot tab.
% * *Transport Analysis*: runs model, prompts for Case description which is added to the listing on the Cases tab and the results can be displayed and plotted from the Q-Plot tab. 
% * *Sediment Properties*: utility to compute a range of sediment properties based in the sediment bulk density and optionally other properties such as temperature and salinity.
% * *Derive Output*: initialises the Derive Output UI to select and define manipulations of the data or call external functions and load the result as new data set.

%% Analysis
% * *Plots*: initialises the Plot UI to select variables and produce various types of plot. The user selects the Case, Dataset and Variable to used, along with the Plot type and any Scaling to be applied from a series of drop down lists, 
% * *Statistics*: initialiss the Statistics UI to select data and run a range of standard statistical methods.

%% Help
% * *Help*: access the online documentation for CoastalTools.

%% See Also
% The <matlab:sdt_open_manual manual> provides further details of setup and 
% configuration of the model.