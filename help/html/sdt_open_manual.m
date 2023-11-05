function sdt_open_manual()
%find the location of the asmita app and open the manual
appinfo = matlab.apputil.getInstalledAppInfo;
idx = find(strcmp({appinfo.name},'SedTools'));
fpath = [appinfo(idx(1)).location,[filesep,'SedTools',filesep,'doc',filesep,'SedTools manual.pdf']];
open(fpath)
