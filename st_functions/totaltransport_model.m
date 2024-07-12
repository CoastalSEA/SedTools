function [results,meta,modeltime] = totaltransport_model(mobj)
%
%-------header-------------------------------------------------------------
% NAME
%   totaltransport_model.m
% PURPOSE
%   run the total load transport model for SedTools App
% USAGE
%    [results,xy,modeltime] = totaltransport_model(mobj) 
% INPUTS
%   mobj - handle to SedTools model instance
% RESULTS
%   results - depth, tidal velocity, bed shear stress, total load transport
%   meta - meta data for additional single valued results
%   modeltime - time of simulation in days
% NOTES
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%
    results =[]; meta = []; modeltime = [];
    %handles to input classes
    tideobj = mobj.Inputs.TideParams;     
    siteobj = mobj.Inputs.SiteParams;  
    rnpobj = mobj.Inputs.RunParams;
    %get input data as structs
    tidespec = getPropertiesStruct(tideobj);
    site = getPropertiesStruct(siteobj); 
    seds = getPropertiesStruct(mobj.Inputs.TransportParams);
    cn = getConstantStruct(mobj.Constants);

    if tideobj.VelocityPhase==0 && siteobj.BedRoughness>0
        %use bed roughness to set phase of velocity (ie k.x+phi) (rads) if
        %not specified as a tidal parameter (i.e. VelocityPhase = 0)
        phi = FrictionVelocityPhase(siteobj,tideobj);    %phase lag of velocity relative to elevation (rads)
        tidespec.VelocityPhase = tidespec.ElevPhase+phi*180/pi; %phase of velocity (ie k.x+phi) (deg)
    end

    %setup time domain
    % tstep = 10;              %timestep (minutes)
    % delt = tstep*60;         %time step(secs) 
    % rundur = 1;              %run duration (days)
    % Tsn  = rundur*24*3600;   %duration of spring-neap cycle (secs)
    % t = 0:delt:Tsn;          %run for spring cycle
    % t = t(1000:ceil(1000+2*12.4*3600/delt));

    %get tidal elevation and velocity
    tide = simple_tide(tidespec,rnpobj.Duration,rnpobj.Tinterval);

    %depth of water over bank
    d = (tide.z-site.zBed);                 
    dmin = 0.1;
    d = d.*(d>dmin);  
    u = tide.u;
    u = u.*(d>0);

    %create arrays based on input wave parameters
    blanks = ones(length(u),1);
    Hs = blanks*site.Hs;
    Tp = blanks*site.Tp; 
    phid = blanks*site.Phid;

    %check for depth limiting condition
    [~,Hs] = hs_break(Hs,Tp,site.BedSlope,d,1,cn.g,1);
    Hrms =  Hs/sqrt(2);   %root mean square wave heights
    
    %temperature/salinity dependent density of water and kinematic viscosity
    [cn.rhow,cn.visc] = fluidprops(seds.Salinity,seds.WaterTemp);
    
    %bed shear stress for combined waves and currents
%     if site.RippleHeight>0
%         %For ripples of height, D and wavelength L, 
%         %an equivalent d50 grainsize of approximately d50 = 12*D^2/L can be used.
%         d90tod50 = seds.D90/seds.D50;
%         seds.D50 = 12*site.RippleHeight^2/site.RippleLength;
%         seds.D90 = seds.D50*d90tod50;
%     end
    d50v = blanks*seds.D50;
    d90v = blanks*seds.D90;
    
    tau = tau_bed(d,d50v,cn.visc,cn.rhow,abs(u),Hrms,Tp,phid);
    taux  = tau.taux; %maximum wave-current induced bed shear stress (N/m2)    
    taux(isnan(taux)) = 0;
    taur = tau.taur;  %root-mean-square wave-current induced bed shear stress (N/m2)
    taur(isnan(taur)) = 0;

    %critical bed shear stress for given sediment size   
    [tau.taucr,~] = tau_crit(seds.BedDensity,seds.D50,seds.PercentMud/100,cn.visc,cn.rhow,cn.rhos);

    qt = svr_transport(d,u,Hs,Tp,d50v,d90v,1/site.BedSlope,cn);
    if isempty(qt), return; end
    netqt = sum(qt,'omitnan');
    [delta,lambda] = sandwave(d,tau.taux,tau.taucr,seds.D50);
    txt1 = sprintf('Erosion threshold = %1.3g N/m^2; Net transport = %1.3g m^3/m/s',tau.taucr,netqt);
    txt2 = sprintf('Wavelength = %1.3g m; Waveheight = %1.3g m',max(lambda),max(delta));
    fprintf('%s\n%s\n',txt1,txt2);

    t = seconds(tide.t);
    modeltime = days(t);
    meta = sprintf('%s; %s',txt1,txt2);
    results = {d,u,taux,taur,qt};
    
    mergestructs = @(x,y) cell2struct([struct2cell(x);struct2cell(y)],[fieldnames(x);fieldnames(y)]);
    data = mergestructs(site,tidespec);
    transport_plot(t,d,u,tau,qt,data);
end

%%
function phi = FrictionVelocityPhase(siteobj,tideobj)
    %use the friction factor based on the tidal current amplitude to set 
    %the velocity phase. Based on the linearisation of the qudratic
    %velocity term viz: 8/(3pi)*f*U/h the friction term for tidal propagation 
    %(see Eq 17 in Friedrichs C T and Aubrey D G, 1994, Journal of 
    %Geophysical Research, 99 (C2), pp. 3321-3336) is:
    ff = 8/3/pi*siteobj.BedRoughness*tideobj.VelocityAmp./(tideobj.MSL0-siteobj.zBed);
    %based on Pillsbury derivation of an ideal estuary, where 
    %cot(phi) = ff/omega
    phi = acot(ff/tideobj.AngularFrequency);  %in rads
end