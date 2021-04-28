function [distribution, classes, stat]=cal_settling(data, temp, varargin)
% Convert the record of the settling columns to grain size
%
%    DISTRIBUTION = CAL_SETTLING(DATA, TEMP, DENSITY, H, RESOLUTION, option)
%    [DISTRIBUTION, CLASSES] = ...
%    [DISTRIBUTION, CLASSES, STAT] = ...
%    [Ws] = CAL_SETTLING([], TEMP, ...)
%
% Inputs  DATA : Filename or vector with raw data from the settling column. To 
%                average several sub-samples, DATA can be a cell array of strings   
%                with the different filenames. In this case TEMP, DENSITY and H may  
%                have the same length than filenames in DATA.
%         TEMP : Temperature (°C)
%         DENSITY : Sediment density (optional, default 2650 kg/m3)
%         H : Height (metres, optional, default 1.768 m)
%         RESOLUTION: Resolution for result, 1 for 1 Phi, 2 for 1/2 Phi, 4 for 1/4  
%                     Phi, 8 for 1/8 Phi, etc. (optional, default 1/8 Phi)
%         option : "nofig" - do not create a figure (to use in loops)
%                  "copyall" - copy to the system clipboard 2-columns matrix  
%                              containing filename, DISTRIBUTION and STAT
%                  "copy" - same as "copyall", but only 2nd column (data column)  
%                  "ixx"  - ignore the xx first seconds (example: "i08" for 8 seconds)
% Output  DISTRIBUTION : WITH 1 output argument, 2-columns matrix with finer limit of 
%                             each grain-size classe and its sediment-mass proportion.
%                        WITH 2 or more output arguments, Sediment mass proportion in 
%                             each class (vertical vector), 
%         CLASSES : 2-columns matrix with lower and upper limits of grain-size in Phi,
%                   the first and last classes includes also all sediments which are 
%                   finer or coarser, respectively.
%         STAT : 2-columns cell matrix with grain-size parameters in Phi,
%                first column: parameter name, second column: value.
%         Ws : If DATA is an empty matrix, only the settling velocities (m/s) used for 
%              the conversion are returned. Format: 2-columns matrix with Ws and
%              corresponding grain-sizes in Phi.

% version 1.17
% Copyright 2005 Urs Neumeier

% The input data DATA must be the reading of the balance (units do not matter)
% recorded with a fixed time-step. The time-step must be written as the constant
% time_step in the beginning of the routine. When the input data are given as in 
% a file, it must be a text file with one weight record per line and nothing else.
%
% The formula used to convert settling velocity to grain-size can be selected 
% be changing the constant after the SWITCH statement in the program.
%
% Often a slight increase in the recorded mass exists in the first seconds, 
% which does not correspond to settling grains but represent some noise. This 
% produces a peak in the coarsest grain-size class. To remove this effect, use 
% the option "ixx". For example, to ignore the 8 first seconds for file data.txt
% with temperature 22:   cal_settling('data.txt',20,'i8')

% Rules followed for the conversion
% 1. The limits of the grain-size classes are set.
% 2. The corresponding settling velocities (cl_ws) are calculated.
% 3. The 7 first records (about 1 seconds) are ignored to avoid perturbation due to 
%    the sediment introduction. The maximum of the 7 first values is considered as 
%    the zero and substracted from all the dataset. 
% 4. If the option ixx is used, then the xx first seconds are ignored, and the first
%    values after that time-interval is considered as zero and substracted from all
%    the dataset.
% 5. Should any part the records be decreasing, then this is corrected.
% 6. The settling velocities are calculated (Ws) for each data point.
% 7. The data are resampled (data2), keeping only the first point, the last point
%    that is zero, and all points that are higher than the previous. Ws is resampled
%    to Ws2 accordingly.
% 8. Applying a linear interpolation on data2 and Ws2, the mass at the classe limits
%    is calculated
% 9. The mass of each class is calculated and the results normalized to 100%.

% To select the statistic parameters printed on the graphic output, change the 
% following variable
% 0 : none,  1 : Folk and Ward (graphical method), 2 : moment method
parameters_on_graphic = 1;

time_step = 0.1591; % time-step of the recorder in seconds
% check input arguments
%if nargin==3 & iscell(varargin{1})
%    varargin=varargin{1}  % recursive call for multiple file
%end
recursive=logical(0);
draw_figure=logical(1);
to_clipboard=0; % 0: do nothing, 1: option copy, 2: option copyall
ixx=0; % time-interval (in seconds) ignored at the beginning, specified with option ixx
nb_data=1; % number of input files (for averaging several sub-samples)
for i=length(varargin):-1:1;
    vari=varargin{i};
    if isequal(vari,'recursive')
        recursive=logical(1);
        varargin(i)=[];
    elseif any(strcmpi(vari,{'nofig','-nofig','nofigure','-nofigure'}))
        draw_figure=logical(0);
        varargin(i)=[];
    elseif any(strcmpi(vari,{'copy','-copy'}))
        to_clipboard=1;
        varargin(i)=[];
    elseif any(strcmpi(vari,{'copyall','-copyall'}))
        to_clipboard=2;
        varargin(i)=[];
    elseif ischar(vari)
        if strncmpi(vari,'i',1)
            ixx=str2num(vari(2:end));
            if isempty(ixx) | ixx<0 | ixx>120
                ixx=0;
            end
        elseif strncmpi(vari,'-i',2)
            ixx=str2num(vari(3:end));
            if isempty(ixx) | ixx<0 | ixx>120
                ixx=0;
            end
        end
        if ixx>0 | strcmpi(vari,'i0')
            if ~recursive, fprintf('The first %g seconds are ignored.\n',ixx); end
            varargin(i)=[];
        end
    end
end
error(nargchk(0,3,length(varargin)))
error(nargchk(2,8,nargin))
if iscellstr(data) & length(data)==1, data=data{1}; end
if iscellstr(data)
    nb_data=length(data);
    for i=1:nb_data
        len(i)=length(data{i});
    end
    v=strvcat(data);
    for i=1:min(len)+2
        if i>min(len) | ~all(v(1,i)==v(:,i))
            filename=data{1}(1:i-1);
            if filename(end)=='_', filename(end)=[]; end
            break
        end
    end
else
    if ischar(data)
    	filename=data;
    	load(data)
    	data=eval(filename(1: min(find(filename=='.'))-1 ));
    else
        filename='';
    end
    if ~isempty(data) & (~isnumeric(data) | min(size(data))~=1 | length(data)<20)
    	error('Incorrect first argument')
    end
end

if nb_data~=1 & prod(size(temp))==1, temp=repmat(temp,1,nb_data); end
if ~isnumeric(temp) | (prod(size(temp))~=nb_data) | all(temp<5) | all(temp> 40)
	error('Incorrect second argument')
end

if length(varargin)>=1 & ~isempty(varargin{1})
    density=varargin{1};
else
	density=2650;
end
if nb_data~=1 & prod(size(density))==1, density=repmat(density,1,nb_data); end
if ~isnumeric(density) | prod(size(density))~=nb_data | all(density<1100) | all(density>22000)
	error('Incorrect third argument')
end

if length(varargin)>=2 & ~isempty(varargin{2})
    h=varargin{2};
else
	h=1.768;
end
if nb_data~=1 & prod(size(h))==1, h=repmat(h,1,nb_data); end
if ~isnumeric(h) | prod(size(h))~=nb_data | all(h<1.0) | all(h>5.0)
	error('Incorrect fourth argument')
end

if length(varargin)>=3 & ~isempty(varargin{3})
    resolution=varargin{3};
	if ~isnumeric(resolution) | prod(size(resolution))~=1 | resolution<1 | resolution>100
		error('Incorrect fifth argument')
	end
else
	resolution=8;
end

if nb_data~=1
    % averaging of several sub-samples
    for i=1:nb_data
        [mdistrib(:,i), classes, method] = ...
            cal_settling(data{i},temp(i),density(i),h(i),resolution,sprintf('i%g',ixx),'recursive');
    end
    distribution=mean(mdistrib')';
else

    % processing of only one sample
    [cl_ws,classes,finest,coarsest,method]=settling_velocity(temp,density,resolution);
    % local function returning the limits of each grain-size classes 
    % and the corresponding settling velocities 

    if isempty(data)				% if data is empty, return Ws as unique output argument
    	distribution=[cl_ws,classes];
    	clear classes
    	return
    end

    % Preprocessing of the data
    if ixx==0
        data=data-max(data(1:7));		% the first 7 data points (~ the first second) are discarded
        data(1:7)=0;
    else
        ixxnb=min(  max( ceil(ixx/time_step) ,2)  ,length(data));
        data=data-data(ixxnb);
        data(1:ixxnb-1)=0;
    end
   
    if any(diff(data)<0)	% force data to be always increasing
    	m=0;
    	for i=1:length(data)
    		if data(i)<m
    			data(i)=m;
    		end
    		m=data(i);
    	end
    end
    
    Ws=h./(([0.1 1:length(data)-1])*time_step);% settling velocity for each data-point (m/s)

    data2=data;	% resample data and Ws
    Ws2=Ws;
    pos=2;
    f=find(data>0);
    if isempty(f)
    	Ws2=Ws(1:2);
    	data2=data(1:2);
    else
    	if f(1)>2
    		data2(pos)=0;
    		Ws2(pos)=Ws(f(1)-1);
    		pos=pos+1;
    	end
    	for i=f(1):length(data)
    		if diff(data(i-1:i)) % true if data(i-1) < data(i)
    			data2(pos)=data(i);
    			Ws2(pos)=Ws(i);
    			pos=pos+1;
    		end
    	end
    	data2(pos)=data2(pos-1);	% add last data point for Ws=0 (t=inf)
    	Ws2(pos)=0;
    	data2(pos+1:end)=[];
    	Ws2(pos+1:end)=[];
    end
    
    if max(data2)==0				% if data is empty or always zero, exit function
    	distribution=[];
    	classes=[];
    	return
    end

    % Conversion of settling-velocities to grain-size
    cl_mass=interp1(Ws2,data2,cl_ws,'linear');	% interpolated data2 for classe limits
    distribution=[cl_mass(1); diff(cl_mass); data2(end)-cl_mass(end)];	% compute mass in each class
    distribution=distribution/data2(end);		% normalize to 100%
    classes=[[classes;finest] [coarsest;classes]];		% complete classes with lower classe-limits
    classes=flipud(classes);			% reverse classes and distribution
    distribution=flipud(distribution);
end

if ~recursive	% do not calculate parameter and draw figure if recursive call
	% calculation of grain-size parameters
	if resolution<16 % parameters are calculated of 1/16th Phi distribution
		[distr_16,cl_16] = cal_settling(data,temp,density,h,16,sprintf('i%g',ixx),'recursive');
		stat=compute_parameters(cl_16,distr_16);
	else
		stat=compute_parameters(classes,distribution);
	end
	if nargout==0
        fprintf('Calculation method : %s\n',method);
		for i=1:size(stat,1)
			fprintf('%7s  %g\n',stat{i,1},stat{i,2});
		end
	end
	% draw figure of the frequency distribution
    if draw_figure
    	plot_grain_size(classes,distribution,stat,parameters_on_graphic,method)
    	if ~isempty(filename)
    		title(filename)
        end
	end
else
    if nargout==3
        stat=method;
    end
end

if to_clipboard ~= 0
    copy_to_clipboard(filename,classes,distribution,stat,to_clipboard)
end


if nargout==1	% delete the unnecessary output arguments
	distribution=[classes(:,1) distribution];
	clear classes
elseif nargout==0
	clear distribution classes
end


function [cl_ws,classes,finest,coarsest,method]=settling_velocity(temp,density,resolution)
% function returning the limits of each grain-size classes and the corresponding settling velocities 
%
% calculation of physical constant for the experiment temperature
% density of distilled water from 5 to 40 °C in kg/m^3
water_densities=[999.97 999.94 999.90 999.85 999.78 999.70 999.61 999.50 ...
		999.38 999.25 999.10 998.94 998.78 998.60 998.41 998.21 997.99...
		997.77 997.54 997.30 997.05 996.79 996.52 996.24 995.95 995.65...
		995.34 995.03 994.71 994.38 994.04 993.69 993.33 992.97 992.60 992.22];
water_density=interp1(5:40,water_densities,temp,'spline');
viscosity_20=1.0020;					% dynamic viscosity in mPa s at 20°C
a= (1.1709*(20-temp) - 0.001827*(temp-20)^2) / (temp + 89.93);
viscosity=10^a * viscosity_20 * 0.001;	% dynamic viscosity in Pa s
g=9.81;									% gravity in m/s^2

% calculation of the limits of grain-size classes
max_grain_size=-2;	
min_grain_size=4.5;	% there will be a class behond max_grain_size and min_grain_size
classes=fliplr(0:-1/resolution:max_grain_size);
classes=[classes(1:end-1) 0:1/resolution:min_grain_size]';	% grain-size classes limits in Phi
cl_m=0.001 * 2.^(-classes);					% grain-size classe limits in metres
finest=classes(end) + 1/resolution;
coarsest=classes(1) - 1/resolution;

% Calculation of the settling-velocity of each classe limit (cl_ws)
% Select one of the following calculation method of the selling velocity 
% be setting the corresponding number after "switch"
%
% Each method produce slightly different results, therefore ALL samples of 
% a study MUST be calculated with the SAME method, so that they can be compared.
% The most recent formula from Soulsby (1997) seems to fit the best with experimental data.
%
switch 4
case 1
	% Gibbs' equation for sphere, the equation uses cgm units
	% R.J. Gibbs, M.D. Matthews and D.A. Link (1971) The relationship between sphere size
	% and settling velocity. Journal of sedimentary Petrology, 41/1, 7-18
    method='Gibbs et al. (1971)';
	g_cgm=g*100;	% m/s^2 -> cm/s^2
	Ps=density/1000;	% kg/m^3 -> g/cm^3
	Pf=water_density/1000;	% kg/m^3 -> g/cm^3
	eta=viscosity*10;	% kg/m/s or Pa s -> poise
	r=cl_m*100/2;		% diameter in m -> radius in cm
	cl_ws = (-3*eta + (9*eta^2 + g_cgm*r.^2*Pf*(Ps-Pf).*(0.015476+0.19841*r)) .^ 0.5) ...
		./ (Pf*(0.011607+0.14881*r));
	cl_ws=cl_ws*0.01;	% cm/s -> m/s
case 2
	% Baba & Komar's modification of Gibbs' equation to get the intermediate diameter of 
	% an elipsoide instead of the diameter of a sphere. For 0.2 < d < 2.0 mm
	% Baba J. & Komar P.D. (1981) Measurements and analysis of settling velocities of 
	% natural quartz sand grains. J. sediment. Petrol., 51, 631-640.
    method='Baba and Komar (1981)';
	g_cgm=g*100;	% m/s^2 -> cm/s^2
	Ps=density/1000;	% kg/m^3 -> g/cm^3
	Pf=water_density/1000;	% kg/m^3 -> g/cm^3
	eta=viscosity*10;	% kg/m/s or Pa s -> poise
	r=cl_m*100/2;		% diameter in m -> radius in cm
	cl_ws = (-3*eta + (9*eta^2 + g_cgm*r.^2*Pf*(Ps-Pf).*(0.015476+0.19841*r)) .^ 0.5) ...
		./ (Pf*(0.011607+0.14881*r));
	cl_ws=0.977*cl_ws.^0.913;	% equation (3) of Baba and Komar
	cl_ws=cl_ws*0.01;	% cm/s -> m/s
	%
case 3
	% van Rijn's equation 3.2.22 (for 100 ~< d ~< 1000 µm) to convert Ws of non-sperical sediment 
	% to sieve diameter
	% L.C. van Rijn (1993) Principles of sediment transport in rivers, estuaries and coastal seas
	% Aqua Publications
    method='van Rijn (1993)';
	viscosity=viscosity/water_density;  % the formula require the kinematic viscosity
	cl_ws = 10*viscosity./cl_m .* ...
		((1 + 0.01*(density/water_density-1)*g*cl_m.^3/viscosity^2) .^0.5 - 1);
case 4
	% Soulsby's equation SC(102) to convert Ws of natural sands to sieve diameter
	% Soulsby R. (1997) Dynamics of marine sands: a manual for practical applications. 
	% Thomas Telford, London, 249 p.
    method='Soulsby (1997)';
	viscosity=viscosity/water_density;  % the formula require the kinematic viscosity
	Dstar=(g*(density/water_density-1)/viscosity^2) ^ (1/3) * cl_m;
	cl_ws = viscosity./cl_m .* ((10.36^2 + 1.049*Dstar.^3).^0.5 - 10.36);
end			






function stat=compute_parameters(cl,distr)
% Compute grain-size parameters
%
% Input 1 argument cl = [finer limit of interval in Phi, sediment mass of interval]
%       2 argument cl = sediment mass of interval (vertical array)
%                  distr = [finer limit, coarser limit] of intervals in Phi
%
% For parameters definition see Folk R.L. & Ward W.C. (1957) 
% Brazos River bar: a study in the significance of grain size parameters. 
% J. sediment. Petrol., 27/1, 3-26.
% version 1.04
if nargin==1
	res=-mean(diff(cl(:,1)));
	distr=cl(:,2);
	cl=[cl(:,1) [cl(2:end,1);cl(end,1)-res]];
end
clm=(mean(cl'))';			% calculation of moment
m1=sum(clm .* distr) / sum(distr);
m2=sum((clm-m1).^2 .* distr) / sum(distr);
m3=sum((clm-m1).^3 .* distr) / sum(distr);
m4=sum((clm-m1).^4 .* distr) / sum(distr);

perc=([5 10 16 25 50 75 84 90 95]./100)';	% calculation of Dx
strPerc={'D5';'D10';'D16';'D25';'D50';'D75';'D84';'D90';'D95'};

Idistr=[0; cumsum(distr)];
Icl=[cl(:,1);cl(end,2)];
working_segment=max(find(Idistr==0)):min(find(Idistr>0.99));
Icl=Icl(working_segment);
Idistr=Idistr(working_segment);
if ~all(diff(Idistr))
	for i=2:length(Idistr)
		if Idistr(i)<=Idistr(i-1)
			Idistr(i)=Idistr(i-1)+5*eps;
		end
	end
	Idistr=Idistr./Idistr(end);
end

Dx=interp1(Idistr,Icl,perc,'linear');

Mz=(Dx(3)+Dx(5)+Dx(7))/3;		% calculation of Folk & Ward's paramenters
s1= -( (Dx(7)-Dx(3))/4 + (Dx(9)-Dx(1))/6.6 );
SK1= -( ...
	(Dx(3)+Dx(7)-2*Dx(5)) / (2*(Dx(7)-Dx(3))) + ...
	(Dx(1)+Dx(9)-2*Dx(5)) / (2*(Dx(9)-Dx(1)))...
	);
Kg=(Dx(9)-Dx(1)) / (2.44*(Dx(6)-Dx(4)));

distr=distr*100/sum(distr);
clay = sum(distr(1:max(find(cl(:,2)>=7.96))));
silt = sum(distr(1:max(find(cl(:,2)>=3.96)))) - clay;
sand = sum(distr(1:max(find(cl(:,2)>=1.01)))) - clay - silt;
coarser = sum(distr) - clay - silt - sand;



stat=[{'Mz',Mz
	'sigma1',s1
	'SK1',SK1
	'Kg',Kg}
	[strPerc num2cell(Dx)]
	{'Moment1',m1
	'Moment2',m2
	'Moment3',m3
	'Moment4',m4
	'clay',clay
	'silt',silt
	'sand',sand
	'coarser',coarser}];



function plot_grain_size(varargin)
%                       (classes,distribution,stat,parameters_on_graphic,method)
% plot the grain-size frequency per weight and the cummulative frequency
% on a figure, and calculate the grain-size parameters.
% version 1.05
if ischar(varargin{end})
    method=varargin{end};
    varargin(end)=[];
end
if isnumeric(varargin{end}) & prod(size(varargin{end}))==1
    parameters_on_graphic=varargin{end};
    varargin(end)=[];
else
    parameters_on_graphic=1; % default value
end
if iscell(varargin{end})
    stat=varargin{end};
    varargin(end)=[];
end
if length(varargin)==1
    classes=varargin{1}(:,1);
	distribution=varargin{1}(:,2);
elseif length(varargin)==2
    classes=varargin{1}(:,1);
	distribution=varargin{2}(:,1);
else 
    error('Incorrect arguments')
end
res=-mean(diff(classes(:,1)));
classes=[classes(:,1) [classes(2:end,1);classes(end,1)-res]];
distribution=distribution/sum(distribution);
dx=abs(diff(classes(1,:)));
figure	% plot the result on a figure
ha1=axes('units','normalized','position',[0.13 0.11 0.75 0.815]);
hb=bar(classes(:,1)-0.5*dx,distribution,1);
set(hb,'facecolor',[0.5 0.5 1.0])
maxXarray=[1 0.5 0.4 0.25 0.2 0.1 0.05 0.025 0.02 0.01];
maxX=maxXarray(end);
for i=2:length(maxXarray)
	if max(distribution)>maxXarray(i)
		maxX=maxXarray(i-1);
		break
	end
end
set(ha1,'Ylim',[0 maxX],'Ytick',0:maxX/10:maxX,'Yticklabel',sprintf('%g%%|',0:maxX*10:maxX*100))
ylabel 'Grain-size frequency'
ha2=axes('position',get(ha1,'position'));
plot([classes(1,1); classes(:,2)], [0; cumsum(distribution)],'r')
set([ha1,ha2],'Xdir','reverse','Xlim',[classes(end,2) classes(1,1)],...
	'Xtick',ceil(classes(end,2)):floor(classes(1,1)))
set(ha2,'YAxisLocation','right','Ylim',[0 1],'Ytick',0:0.1:1,...
	'YtickLabel',sprintf('%3g%%|',0:10:100),'color','none')
ylabel 'Cummulative frequency distribution'
xlabel 'Grain size (Phi)'
grid on
if exist('method','var')
    Xlim=get(gca,'xlim');
    Ylim=get(gca,'ylim');
    ht2=text(Xlim(2)+diff(Xlim)*0.12,Ylim(1)-diff(Ylim)*0.09,...
        ['Calculation method: ',method],'FontUnits','points','FontSize',7);
end
if  exist('stat','var') & iscellstr(stat(:,1)) 
    if exist('parameters_on_graphic','var') & parameters_on_graphic==0
        % do nothing    
    elseif exist('parameters_on_graphic','var') & parameters_on_graphic==1
        if any(strcmpi('D50',stat(:,1))) & any(strcmpi('Mz',stat(:,1))) & ... 
                any(strcmpi('sigma1',stat(:,1)))
        	D50=stat{strcmpi('D50',stat(:,1)),2};
        	Mz=stat{strcmpi('Mz',stat(:,1)),2};
        	s1=stat{strcmpi('sigma1',stat(:,1)),2};
        	ch=sprintf('D_5_0  %4.3f\nMz  %4.3f\n\\sigma_1  %4.3f',D50,Mz,s1);
        	if any(strcmpi('SK1',stat(:,1))) & any(strcmpi('Kg',stat(:,1)))
        		SK1=stat{strcmpi('SK1',stat(:,1)),2};
        		Kg=stat{strcmpi('Kg',stat(:,1)),2};
        		ch=sprintf('%s\nSK_1%4.3f\nK_G  %4.3f',ch,SK1,Kg);
        	end
        end
    elseif exist('parameters_on_graphic','var') & parameters_on_graphic==2
        if any(strcmpi('Moment1',stat(:,1))) & any(strcmpi('Moment2',stat(:,1))) & ... 
                any(strcmpi('Moment3',stat(:,1))) & any(strcmpi('Moment4',stat(:,1))) ...
                & any(strcmpi('D50',stat(:,1))) 
        	D50=stat{strcmpi('D50',stat(:,1)),2};
        	M1=stat{strcmpi('Moment1',stat(:,1)),2};
        	M2=stat{strcmpi('Moment2',stat(:,1)),2};
        	M3=stat{strcmpi('Moment3',stat(:,1)),2};
        	M4=stat{strcmpi('Moment4',stat(:,1)),2};
        	ch=sprintf(['D_5_0  %4.3f\nx  %4.3f\n\\sigma^2 %4.3f\n',...
                    '\\alpha^3 %4.3f\n\\kappa^4  %4.3f'],D50,M1,M2,M3,M4);
        end
    end
    if exist('ch','var') & ~isempty(ch)    
    	hr=rectangle('position',[classes(1,1)-0.15,0.97,0.01,0.01],...
    		'FaceColor',[0.99 0.99 0.99],'EdgeColor','none','tag','cal_settling rectangle');
    	ht=text(classes(1,1)-0.15,0.97,ch,'VerticalAlignment','top',...
    		'tag','cal_settling text');
    	text_extend=get(ht,'extent');
    	text_extend(1)=text_extend(1)-text_extend(3);
    	set(hr,'position',text_extend);
    	callback_str=['ex=get(findobj(gcbo,''tag'',''cal_settling text''),''extent'');',...
    		'ex(1)=ex(1)-ex(3);',...
    		'set(findobj(gcbo,''tag'',''cal_settling rectangle''),''position'',ex)'];
    	set(gcf,'ResizeFcn',callback_str)
    end
end



function copy_to_clipboard(filename,classes,distribution,stat,option)
if option==1 %copy
    ch=[filename sprintf('\n%18.10f',distribution) sprintf('\n%18.10f',stat{:,2})];
else %copy all
    stat1=stat';
    ch=[sprintf('coarser than\t%s',filename) ...
            sprintf('\n%10.5f\t%18.10f',[classes(:,1) distribution]') ...
            sprintf('\n%s\t%18.10f',stat1{:})];
end
CLIPBOARD('copy',ch)