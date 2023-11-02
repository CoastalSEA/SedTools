function transport_plot(t,d,u,tau,qt,meta)
%
%-------header-------------------------------------------------------------
% NAME
%   transport_plot.m
% PURPOSE
%   plot results from the total transport model in the SedTools App
% USAGE
%       transport_plot(t,d,u,tau,qt)
% INPUTS
%   t - model time (s)
%   d - water depth at site (m)
%   u - tidal velocity at site (m/s)
%   tau - struct of wave-current bed shear stresses using:
%       taum - mean wave-current induced bed shear stress (N/m2)
%       taux - maximum wave-current induced bed shear stress (N/m2)
%       taur - root-mean-square wave-current induced bed shear stress (N/m2)
%   qt - total load transport (m2/s)
%   meta - struct of run input data (optional)
% RESULTS
%   plot of results
%
% Author: Ian Townend
% CoastalSEA (c) Nov 2023
%--------------------------------------------------------------------------
%    
    xlimit = 700;
    time = hours(t);
    netqt = sum(qt,'omitnan');

    figure('Tag','PlotFig',...
        'Units','normalized','Position',[0.2 0.4 0.38 0.34], ...
        'Resize','on','HandleVisibility','on');
    %upper plot---------------------------------------------------------
    subplot(2,1,1);
    yyaxis left
    plot(hours(t),d);
    %xlabel('Time (hours)')
    if time>500
        xlim([0,xlimit]);
    end
    ylabel('Water depth (m)')
    yyaxis right
    plot(time,u);
    ylabel('Velocity (m/s)')

    %lower plot----------------------------------------------------------
    subplot(2,1,2);
    yyaxis left
    plot(time,tau.taux,'DisplayName','Max shear stress');
    hold on   
    plot(time,tau.taur,'--g','DisplayName','RMS shear stress');
    plot(time,tau.taum,':m','DisplayName','Mean shear stress');
    plot([time(1),time(end)],[tau.taucr,tau.taucr],'-.r','DisplayName','Crit. shear stress');
    hold off
    xlabel('Time (hours)')
    if time>500
        xlim([0,xlimit]);
    end
    ylabel('Bed shear stress (N/m2)')
    yyaxis right
    plot(time,qt,'DisplayName','Total load transport');
    ylabel('Total load transport (m2/s)')
    legend
    %added text----------------------------------------------------------
    txt = sprintf('Net Transport over period %1.3g (m3/m-width), Critical shear stress %1.3g (Pa)',...
                                                         netqt,tau.taucr);
    uicontrol(...
    'Style','text',...
    'String', txt,...
    'HorizontalAlignment', 'left',...
    'Units','normalized', ...
    'Position', [0.15 0.45 0.7 0.03],...
    'Tag','Xtext');
    if nargin<6
        sgtitle('Bed stress and Total load transport')
    else
        titletxt = sprintf('Wave height=%.1f, Wave period=%.1f\nTidal amp=%.1f, Velocity phase=%.1f',...
                           meta.Hs,meta.Tp,meta.TidalAmp,meta.VelocityPhase);
        sgtitle(titletxt,'FontSize',12)
    end
end