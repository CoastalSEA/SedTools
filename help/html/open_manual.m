function open_manual()
%find the location of the asmita app and open the manual
% appinfo = matlab.apputil.getInstalledAppInfo;
toolboxes = matlab.addons.toolbox.installedToolboxes;
% idx = find(strcmp({toolboxes.Name},'ModelUI'));
idx = find(strcmp({appinfo.name},'SedTools'));
fpath = [toolboxes(idx(1)).location,'/doc/SedTools manual.pdf'];
open(fpath)
