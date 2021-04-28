function [grainsize,results,stats,method] = settling_column(data,inp)
% NAME
%   settling_column.m
% PURPOSE
%   Convert the record of the settling columns to grain size
% USAGE
%   [results,stats] = settling_column(data,inp,filename)
% INPUTS
%   data - vector of raw data from the settling column. To average several 
%          sub-samples, data can have multiple columns. For multiple data
%          sets temp, density and column height can be single values or an
%          array of the same size
%   inp - structure containing:
%           temp - temperature (°C)
%           density - Sediment density (optional, default 2650 kg/m3)
%           h - Height (metres, optional, default 1.768 m)
%           resolution - Resolution for result, 1 for 1 Phi, 2 for 1/2 Phi, 4 for 1/4  
%                     Phi, 8 for 1/8 Phi, etc. (optional, default 1/8 Phi)
%           ixx - ignore the xx first seconds (example: "i08" for 8 seconds)
% RESULTS
%   results - 2-column cell matrix of grain-size analysis output
%                first column: central value of each class, 
%                second column: sediment mass proportion in each class                             
%   stats - 2-column cell matrix for grain-size parameters,
%                first column: parameter name, second column: value.
% NOTES
%  based on cal_settling by Urs Neumeier, 2005, http://neumeier.perso.ch/matlab/cal_settling.html
% AUTHOR
%   Ian Townend
% COPYRIGHT
%   CoastalSEA (c) Feb 2019
%
    temp = inp.Temperature;
    density = inp.SedimentDensity;
    h = inp.ColumnHeight;
    ixx = inp.OmitPeriod;
    resolution = inp.Resolution;
    settling = inp.SettlingMethod;
    %parameters_on_graphic = inp.PlotStatOption;
    time_step = inp.RecTimeStep; % time-step of the recorder in seconds
    units = inp.OutputUnits; 

    nrec = size(data,2);
    if nrec>1  %force other variables to the same dimension as nrec
        temp = checklength(temp,nrec);
        density = checklength(density,nrec);
        h = checklength(h,nrec);
    end
    %get the distribution for each input data set    
    for i=1:nrec
        [distribution(:,i),classes,method] = ...
            cal_settling(data{i},temp(i),density(i),h(i),resolution,...
                                                ixx,settling,time_step);
    end
    %take mean if data input has multiple data sets
    if nrec >1
        distribution= mean(distribution,2);
    end

    stats = compute_parameters(classes,distribution);
    %adjust units if required
    grainsize = [classes(:,1);classes(end,2)];
    grainsize = classes(:,1)+diff(grainsize)/2;
    if strcmp(units,'mm')
        grainsize = 2.^(-grainsize);
        stats(1:2,2) = num2cell(2.^(-cell2mat(stats(1:2,2))));
        stats(5:14,2) = num2cell(2.^(-cell2mat(stats(5:14,2))));
    end
    grainsize = {grainsize(distribution>0)};
    sedfreq = distribution(distribution>0)*100;
    cumfreq = cumsum(sedfreq);
    results = {cumfreq,sedfreq}; %convert distribution to %
end
%%
function [distribution,classes,method] = cal_settling(data,temp,density,...
                                       h,resolution,ixx,settling,time_step)
    %??
    [cl_ws,classes,finest,coarsest,method]=settling_velocity(temp,density,resolution,settling);
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
    cl_mass=interp1(Ws2,data2,cl_ws,'linear');	   %interpolated data2 for classe limits
    distribution=[cl_mass(1); diff(cl_mass); data2(end)-cl_mass(end)];	%compute mass in each class
    distribution=distribution/data2(end);		   %normalize to 100%
    classes=[[classes;finest] [coarsest;classes]]; %complete classes with lower classe-limits
    classes=flipud(classes);			           %reverse classes and distribution
    distribution=flipud(distribution);
end
%%
function [cl_ws,classes,finest,coarsest,method]=settling_velocity(temp,density,resolution,settling)
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
    switch settling
        case 1
            % Gibbs' equation for sphere, the equation uses cgm units
            % R.J. Gibbs, M.D. Matthews and D.A. Link (1971) The relationship between sphere size
            % and settling velocity. Journal of sedimentary Petrology, 41/1, 7-18
            method={'Gibbs et al. (1971)'};
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
            method={'Baba and Komar (1981)'};
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
            method={'van Rijn (1993)'};
            viscosity=viscosity/water_density;  % the formula require the kinematic viscosity
            cl_ws = 10*viscosity./cl_m .* ...
                ((1 + 0.01*(density/water_density-1)*g*cl_m.^3/viscosity^2) .^0.5 - 1);
        case 4
            % Soulsby's equation SC(102) to convert Ws of natural sands to sieve diameter
            % Soulsby R. (1997) Dynamics of marine sands: a manual for practical applications.
            % Thomas Telford, London, 249 p.
            method={'Soulsby (1997)'};
            viscosity=viscosity/water_density;  % the formula require the kinematic viscosity
            Dstar=(g*(density/water_density-1)/viscosity^2) ^ (1/3) * cl_m;
            cl_ws = viscosity./cl_m .* ((10.36^2 + 1.049*Dstar.^3).^0.5 - 10.36);
    end
end
%%
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
    working_segment=find(Idistr==0, 1, 'last' ):find(Idistr>0.99, 1 );
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
    clay = sum(distr(1:find(cl(:,2)>=7.96, 1, 'last' )));
    silt = sum(distr(1:find(cl(:,2)>=3.96, 1, 'last' ))) - clay;
    sand = sum(distr(1:find(cl(:,2)>=1.01, 1, 'last' ))) - clay - silt;
    coarser = sum(distr) - clay - silt - sand;
    %
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
end
%%
function var = checklength(var,nrec)
        if length(var)~=nrec
            var = repmat(var,1,nrec);
        end
end
