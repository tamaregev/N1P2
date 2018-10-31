function [peak_smpls,  peak_amps, peak_times, manuals ] = PeakDetection_N1P2(grandMatrix, whichSubjects, t, electrodeName, pwin, positive_peak, select_largest, plotflag, block )
%Oct 22 2018 - Tamar
% adapted from PeakDetection in MMNchroma
    LP = false;
    if strcmp(electrodeName,'GFP')
        electrode = 'GFP';
    else
        if strcmp(electrodeName,'central cluster')
            electrodeNames = {'FC1','FCz','FC2',...
                            'C1','Cz','C2',...
                            'CP1','CPz','CP2'};
                                load('L:\Experiments\MMNchroma\Analysis\Properties')
            [ electrode, ~ ] = ChannelName2Number( Properties, electrodeNames );
        else
            load('L:\Experiments\MMNchroma\Analysis\Properties')
            [ electrode, ~ ] = ChannelName2Number( Properties, electrodeName );
        end
    end
    bslwin = [-0.1 0];
    srate = 512;
    range = nan(1,2);[~, range(1)] = min(abs(t-pwin(1)));[~, range(2)] = min(abs(t-pwin(2)));
    [~, onsetsamp] = min(abs(t-bslwin(1)));
    peak_smpls = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4));
    peak_times = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4));
    peak_amps = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4));
    manuals = nan(size(peak_amps));
    failed_sph = cell(size(peak_amps));
    failed_sps = nan(size(peak_amps));
    
    for s = whichSubjects
        h=ERPfigure;
        %set(h,'Position',[80,80,900,600]);
        set(h,'units','normalized','outerposition',[0 0 1 1])
        
        sbplots = randperm(20);%in order to be blind to the condition.
        isp = 0;
        for con = 1:size(grandMatrix,3)
            for prevcon = 1:size(grandMatrix,4)
                
                %take care of electrode
                if isstring(electrode)
                    if strcmp(electrode,'GFP')
                        %calculate GFP
                        %wave = ...
                    end
                elseif numel(electrode)==1
                    wave = grandMatrix(:,s,con,prevcon,electrode);
                elseif numel(electrode)>1
                    wave = mean(grandMatrix(:,s,con,prevcon,electrode),5);
                end
                if LP 
                    wave=LPF(wave,srate,20);
                end
                if ~isnan(wave)
                    if ~(all(~wave))
                        isp = isp+1;
                        sph = subplot(4,5,sbplots(isp));
                        disp(['Elect = ' electrodeName ', subj = ' num2str(s) ', con=' num2str(con) ', prevcon=' num2str(prevcon)])
                        %peak selection
                        select_manual = false;
                        [peak_smpl, manual] = PeakSelect_N1P2(wave, t, [range(1), range(end)], positive_peak, num2str(s), select_manual, select_largest, plotflag, onsetsamp);            
                        
                        if sbplots(isp)==20, xlabel('ms'); end
                        
                        peak_smpls(s,con,prevcon) = peak_smpl;

                        manuals(s,con,prevcon) = manual;
                        if isnan(manual)
                            failed_sph{s,con,prevcon} = sph;
                            failed_sps(s,con,prevcon) = 1;
                        end
                        if ~isnan(peak_smpl)
                            peak_amps(s,con,prevcon) = wave(peak_smpl);
                            peak_times(s,con,prevcon) = t(peak_smpl);
                        else
                            peak_amps(s,con,prevcon) = nan;
                            peak_times(s,con,prevcon) = nan;
                        end
                    end
                end
            end
        end
        if positive_peak
            whichpeak = 'P2';
        else
            whichpeak = 'N1';
        end
        suptitle(['Subject ' num2str(s) ', ' whichpeak ', Block ' block])
        pause
        if ~all(all(isnan(squeeze(failed_sps(s,:,:)))))
            %manually select one by one
            [cons,prevcons]=find(~isnan(squeeze(failed_sps(s,:,:))));
            select_manual = true;
            
            for ii=1:length(cons)
                con = cons(ii);
                prevcon = prevcons(ii);
                disp(['Elect = ' electrodeName ', subj = ' num2str(s) ', con=' num2str(con) ', prevcon=' num2str(prevcon)])
                %take care of electrode
                if isstring(electrode)
                    if strcmp(electrode,'GFP')
                        %calculate GFP
                        %wave = ...
                    end
                elseif numel(electrode)==1
                    wave = grandMatrix(:,s,con,prevcon,electrode);
                elseif numel(electrode)>1
                    wave = mean(grandMatrix(:,s,con,prevcon,electrode),5);
                end
                if LP 
                    wave=LPF(wave,srate,20);
                end
                subplot(failed_sph{s,con,prevcon})
                [peak_smpl, manual] = PeakSelect_N1P2(wave, t, [range(1), range(end)], positive_peak, num2str(s), select_manual, select_largest, plotflag, onsetsamp);            
                if sbplots(isp)==20, xlabel('ms'); end

                peak_smpls(s,con,prevcon) = peak_smpl;
                manuals(s,con,prevcon) = manual;
                if ~isnan(peak_smpl)
                    peak_amps(s,con,prevcon) = wave(peak_smpl);
                    peak_times(s,con,prevcon) = t(peak_smpl);
                else
                    peak_amps(s,con,prevcon) = nan;
                    peak_times(s,con,prevcon) = nan;
                end
            end
            pause
        end
    end

end

