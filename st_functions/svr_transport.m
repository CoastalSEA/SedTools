function qt = svr_transport(d,U,Hs,Tp,d50,d90,beta,cn) 
%
%-------header-------------------------------------------------------------
% NAME
%   svr_transport.m
% PURPOSE
%   function to calculate the total load transport using the equations of
%   Soulsby-van Rijn (p183 Eq(136) in Dynamics of Marine Sands)
% USAGE
%    qt = svr_transport(d,U,Hs,Tp,d50,d90,beta,cn) 
% INPUTS
%   d - water depth (m)
%   U - depth averaged current velocity (m/s)
%   Hs - significant wave height (m)
%   Tp - peak wave period (s)
%   d50 - median sediment grain size (m)
%   d90 - 90 percentile sediment grain size (m)
%   beta - slope of bed in streamwise direction (+ve if flow is uphill) 
%   cn - constants struct (g, rhow, rhos, visc)
% RESULTS
%   qt - total load sediment transport 
% NOTES
%   develped for use on Flora Bank project as totaltransport.m
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%
    depthmask = (d==0);
    %critical threshold velocity
    if length(d50)>1
        id1 = find(d50>=0.0001 & d50<=0.0005);
        id2 = find(d50>0.0005 & d50<=0.002);
        Ucr(id1,1) = 0.19*d50(id1).^0.1.*log10(4*d(id1)./d90(id1));
        Ucr(id2,1) = 8.5*d50(id2).^0.6.*log10(4*d(id2)./d90(id2));
        Ucr(depthmask) = 0;
    else
        if d50>=0.0001 && d50<=0.0005
            Ucr = 0.19*d50.^0.1.*log10(4*d./d90);
        elseif d50>0.0005 && d50<=0.002
            Ucr = 8.5*d50.^0.6.*log10(4*d./d90);
        else
            warndlg('d50 outside valid range in svr_transport');
            qt = [];
            return;
        end
    end
    
    %root mean square wave orbital velocity
    Urms = rmswaveorbitalvelocity(Hs,Tp,d,cn.g);
    Urms(depthmask) = 0;
    
    %specific density of sediment
    s = cn.rhos/cn.rhow;
    %bed roughness length for rippled beds defined in method as 6mm
    z0 = 0.006;
    %dimensionless grain size
    Dst   = d50*(cn.g*(s-1)/cn.visc^2)^(1/3);
    %drag coefficient due to current alone
    CD = (0.4./(log(d/z0)-1)).^2;
    %bed load transprot coefficient
    Asb = 0.005*d.*(d50./d).^(1.2)./((s-1)*cn.g*d50).^1.2;
    Asb(depthmask) = 0;
    %suspended load transport coefficient
    Ass = 0.012*d50.*Dst.^-0.6./((s-1)*cn.g*d50).^1.2;
    %total coefficient
    As = Asb+Ass;
    
    %total transport per unit width of bed (m3/m/s)
    slopefactor = 1-1.6*tan(beta);
    threshold = sqrt(U.^2+0.018./CD.*Urms.^2)-Ucr;
    threshold(depthmask) = 0;
    thresholdfactor = (threshold.*(threshold>0)).^2.4;
    qt = As.*U.*thresholdfactor.*slopefactor; 
end

