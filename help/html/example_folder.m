function example_folder()
%find the location of the example folder and open it
%toolbox
% idx = find(strcmp({toolboxes.name},'muitoolbox'));
% fpath = [toolboxes(idx(1)).location,'/example'];
%app
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'SedTools'));
fpath = [appinfo(idx(1)).location,'/example'];

winopen(fpath)