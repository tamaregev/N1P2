function [peak_smpls, manuals] = PeakSelect_N1P2_3(wave, t, ranges, positives, sName, select_manual, select_largest, plotflag, onsetsamp, bold)
% Tamar Oct 22 - adapted from PeakSelect
% - added t for plotting
% - added select_manual true or false to activate manual select or not. 

% Peak Select is a GUI that allows the user to select manually a local maxima or
% minima within a given range, in case that the algorithm fails. If there is more than one peak 
% within range, the GUI will pop up - a plot and a selection box will be
% displayed and the user will be able to select the peak by its number.
% if there is no peak - the user can select manually
% 
% INPUT:
%
% wave: a vector of voltage values 
%
% range: range in samples (from 1 to n samples in wave) to search for peaks
%
% positive: a boolean with true for positive peaks or false for negative
% peaks
%
% sName: subject's number for display in title of plot
%
% largest:  a boolean with true for picking largest peak (absolute value)
% if multiple
%
% plotflag: true if wanna plot when only 1 peak found, false if not.
%
% onsetsamp: for plotting the time of onset as well. his helps for the
% manual selection
%
% onsetsamp:
%
% OUTPUT:
%
% peak_smpl: the peak sample selected by user
%
% Author: Doron Ariav (email: ariavdoron@gmail.com)

% 8/3/2016 - Tamar: changed window style of inputdlg to normal in order to
% be able to select manually the peaks in case there are none
% Jan 4 2017 Tamar: added inputs: select_largest, plotflag, onsetsamp
% Apr 26 2017 Tamar: replaced close all with close
srate = 512;

if nargin < 6
    select_manual = true;
    select_largest = false;
    plotflag = true;
    onsetsamp = 0;
    bold = false;
end

if nargin < 7
    select_largest = false;
    plotflag = true;
    onsetsamp = 0;
    bold = false;
end

if nargin < 8
    plotflag = true;
    onsetsamp = 0;
    bold = false;
end

if nargin < 9
    onsetsamp = 0;
    bold = false;
end

if nargin < 10
    bold = false;
end

peak_smpls = nan(size(ranges));
manuals = nan(size(ranges));
Colors = {'r','g'};

%loop for N1 and P2
for ipeak=1:length(positives)

    positive=positives(ipeak);
    if ~positive
        wave = wave*-1;
    end

    range = ranges{ipeak};
    
    [pks, locs] = findpeaks(wave);
    local_pks = (locs >= range(1) & locs <= range(2));
    n_peaks = sum(local_pks);
    local_pks = find(local_pks);

    if select_largest
        if n_peaks>1
            %remain with larger peak
            [~, I]=max(wave(locs(local_pks)));
            local_pks=local_pks(I);
            n_peaks=1;
        end
    end

    txt = 'P2';
    if ~positive
        wave = wave*-1;
        pks = pks*-1;
        txt = 'N1';
    end

    options.WindowStyle = 'normal';
    manual = 0;
    if n_peaks < 1
        plot(t,wave);
        xlim([-100 400])
        hold on;
        y = ylim;
        switch ipeak
            case 1 %N1
                rectangle('Position',[t(range(1)),y(1), t(range(2))-t(range(1)), y(2)-y(1)],'FaceColor',[0.8 0.8 1])
            case 2 %P2
                 rectangle('Position',[t(range(1)),y(1), t(range(2))-t(range(1)), y(2)-y(1)],'FaceColor',[1 0.8 0.8])
        end
        if bold
            plot(t,wave,'linewidth',3,'color','k');
        else
            plot(t,wave,'color','b');
        end
        if onsetsamp 
            plot([t(onsetsamp) t(onsetsamp)], [y(1) y(2)], '-.k');
        end
        plot([-100 400],[0 0],'-.k')
        if select_manual
            selection = {};
            while isempty(selection)
                [x1, y1]=ginput(1);
                x1=round(x1);
                plot(x1,y1,'x','MarkerSize',10,'Color','r');
                selection = inputdlg('Press OK to approve selection','manual selection',1,{num2str(x1)});     
            end
            [~, peak_smpl] = min(abs(t-x1));
            manual = manual + 1;disp('manual')
        else
            manual = nan;
            peak_smpl = nan;
        end

    elseif n_peaks > 1 %if select_largest is true than this is never the case
        plot(wave);
        xticklabels({'-100','0','100','200','300','400','500','600'})
        xlabel('ms')
        hold on;
        y = ylim;
        rectangle('Position',[range(1),y(1), range(2)-range(1), y(2)-y(1)],'FaceColor',[0.5 0.5 0.5])
        plot(locs(local_pks),pks(local_pks),'x','MarkerSize',10,'Color','r');
        for i = 1:length(local_pks)
            text(locs(local_pks(i)),pks(local_pks(i))+0.25,num2str(i));
        end
        title(strcat(sName,' : ',txt));
        pos = get(gcf,'position');
        set(gcf,'position', [1 pos(2), pos(3), pos(4)]);
        user_selection = inputdlg('select peak number','',1,{'1'}); 
        close;
        user_selection = str2double(user_selection); 
        peak_smpl = locs(local_pks(user_selection));
    else % one peak found in range
        peak_smpl = locs(local_pks);
        if plotflag
            hold on;
            switch ipeak
                case 1 %N1
                    plot(t,wave);
                    y = ylim;
                    rectangle('Position',[t(range(1)),y(1), t(range(2))-t(range(1)), y(2)-y(1)],'EdgeColor',Colors{ipeak})
                    if bold
                        plot(t,wave,'linewidth',3,'color','k');
                    else
                        plot(t,wave);
                    end
                    xlim([-100 400])
                    if onsetsamp 
                        plot([0 0], [y(1) y(2)], '-.k');
                    end
                    plot([-100 400],[0 0],'-.k')
                    plot(t(locs(local_pks)),pks(local_pks),'x','MarkerSize',10,'Color',Colors{ipeak});
                    for i = 1:length(local_pks)
                        text(locs(local_pks(i)),pks(local_pks(i))+0.5,txt,'Color',Colors{ipeak});
                    end
                case 2 %P2
                    rectangle('Position',[t(range(1)),y(1), t(range(2))-t(range(1)), y(2)-y(1)],'EdgeColor',Colors{ipeak})
                    plot(t(locs(local_pks)),pks(local_pks),'x','MarkerSize',10,'Color',Colors{ipeak});
                    for i = 1:length(local_pks)
                        text(locs(local_pks(i)),pks(local_pks(i))+0.5,txt,'Color',Colors{ipeak});
                    end
            end
        end
    end
   peak_smpls(ipeak) = peak_smpl;
   manuals(ipeak) = manual;
end
%close;
end

