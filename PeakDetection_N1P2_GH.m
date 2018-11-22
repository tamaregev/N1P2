function [peak_smpls,  peak_amps, peak_times, grandcon_as_median, grandcon_peak_times, grandcon_peak_amps] = PeakDetection_N1P2_GH(grandMatrix, whichSubjects, t, electrodeName, dt, pwins, positive_peaks, select_largest, plotflag, pauseflag, block )
%Oct 22 2018 - Tamar
% adapted from PeakDetection_N1P2_2 in N1P2
%TODO grandcon_peak_times is redundant with peak_times
%Nov 14 2018 - Tamar - changed name from PeakDetection_N1P2_3 to PeakDetection_N1P2_GH
%Nov 14 2018 - Tamar - added checking whether a specific subject
%participated in this block - in MMNchroma few subjects didn't do block
%classicctrl

Colors = {'r','b'};
for ipeak = 1:length(pwins)
    pwin = pwins{ipeak};
    range = nan(1,2);
    [~, range(1)] = min(abs(t-pwin(1)));
    [~, range(2)] = min(abs(t-pwin(2)));
    ranges{ipeak}=range;
end
clear range

srate = 512;
ds = round(srate*dt/1000); %calculate window in samples from time    

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
[~, onsetsamp] = min(abs(t-bslwin(1)));
peak_smpls = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4),2);
peak_times = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4),2);
peak_amps = nan(size(grandMatrix,2),size(grandMatrix,3),size(grandMatrix,4),2);
%if we want to save the grand peaks:
grandcon_peak_smpls = nan(size(grandMatrix,2),size(grandMatrix,3),2);
grandcon_peak_amps = nan(size(grandMatrix,2),size(grandMatrix,3),2);
grandcon_peak_times = nan(size(grandMatrix,2),size(grandMatrix,3),2);
grandcon_as_median = zeros(size(grandMatrix,2),size(grandMatrix,3),2);
%     
    manualss = nan(size(peak_amps,1),size(peak_amps,2),2);
    failed_sph = cell(size(peak_amps,1),size(peak_amps,2),2);
    failed_sps = nan(size(peak_amps,1),size(peak_amps,2),2);
    
    for s = whichSubjects
        %check whether not all grandcon is NaN, meaning that this subject
        %didn't do this condition at all
        sGrand = grandMatrix(:,s,:,:,electrode);
        if ~all(all(all(isnan(sGrand))))
            h=ERPfigure;
            %set(h,'Position',[80,80,900,600]);
            set(h,'units','normalized','outerposition',[0 0 1 1])

            conord = randperm(size(grandMatrix,3));%randomly order conditions
            isp=0;
            for con = 1:size(grandMatrix,3)
                isp=isp+1;
                sph = subplot(5,5,conord(isp));
                %calculate wave of all con:
                if isstring(electrode)
                    if strcmp(electrode,'GFP')
                        %calculate GFP
                        %wave = ...
                    end
                elseif numel(electrode)==1
                    wave = nanmean(grandMatrix(:,s,con,:,electrode),4);
                elseif numel(electrode)>1
                    wave = nanmean(mean(grandMatrix(:,s,con,:,electrode),5),4);
                end
                %peak selection grand con  
                select_manual = false;
                bold = true;
                [peak_smpl2, manual2] = PeakSelect_N1P2_3(wave, t, ranges, positive_peaks, num2str(s), select_manual, select_largest, plotflag, onsetsamp, bold);                        

                grandcon_peak_smpls(s,con,:) = peak_smpl2;
                manualss(s,con,:) = manual2;%manual is nan if failed, 1 if manually selected and 0 if all good automatically
                %TODO - decide whether to keep manualss variable or replace
                %with grandcon_as_median
                for ipeak = 1:length(positive_peaks)%check for N1 and P2
                    if isnan(manual2(ipeak))
                        failed_sph{s,con,ipeak} = sph;
                        failed_sps(s,con,ipeak) = 1;
                    end
                    if ~isnan(peak_smpl2(ipeak))
                        grandcon_peak_amps(s,con,ipeak) = wave(peak_smpl2(ipeak));
                        grandcon_peak_times(s,con,ipeak) = t(peak_smpl2(ipeak));
                    else
                        grandcon_peak_amps(s,con,ipeak) = nan;
                        grandcon_peak_times(s,con,ipeak) = nan;
                    end
                end

                %calc and plot prevconds - 
                if all(isnan(failed_sps(s,con,:)))
                %same as: if ~any(isnan(manualss(s,con,:))
                    sbplots = size(grandMatrix,3)*(1:(size(grandMatrix,4)-1))+conord(isp);
                    sbplots = sbplots(randperm(length(sbplots)));
                    ispp=0;
                    for prevcon = 1:size(grandMatrix,4)
                        %calculate 'wave':
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
                        %calc peak and plot:
                        if ~isnan(wave)
                            if ~(all(~wave))
                                ispp = ispp+1;
                                spph = subplot(5,5,sbplots(ispp));
                                disp(['Elect = ' electrodeName ', subj = ' num2str(s) ', con=' num2str(con) ', prevcon=' num2str(prevcon)])
                                %calc and plot N1-P2 amp
                                for ipeak = 1:length(positive_peaks)%run on N1 and P2
                                    smpl = grandcon_peak_smpls(s,con,ipeak);
                                    peak_amp = mean(wave((smpl-ds):(smpl+ds)));                        
                                    %plot wave
                                    hold on;
                                    plot(t,wave);
                                    y = ylim;
                                    rectangle('Position',[t(smpl-ds),y(1), t(smpl+ds)-t(smpl-ds), y(2)-y(1)],'FaceColor',[0.9 0.9 0.9])
                                    plot(t,wave);
                                    xlim([-100 400])
                                    plot([0 0], [y(1) y(2)], '-.k');
                                    plot([-100 400],[0 0], '-.k')
                                    plot([t(smpl-ds) t(smpl+ds)],[peak_amp peak_amp],'Color',Colors{ipeak},'linewidth',2)
                                    peak_smpls(s,con,prevcon,ipeak) = smpl;
                                    if ~isnan(smpl)
                                        peak_amps(s,con,prevcon,ipeak) = peak_amp;
                                        peak_times(s,con,prevcon,ipeak) = t(smpl);
                                    else
                                        peak_amps(s,con,prevcon,ipeak) = nan;
                                        peak_times(s,con,prevcon,ipeak) = nan;
                                    end                            
                                end

                                if sbplots(ispp)==25, xlabel('ms'); end

                            end
                        end
                    end
                end
            end
            if ~all(all(isnan(failed_sps(s,:,:))))
                %compute one by one ; as median of other conditions
                 [cons,ipeaks]=find(~isnan(squeeze(failed_sps(s,:,:))));

                 for ii=1:length(cons)
                     con = cons(ii);
                     for ip = 1:length(ipeaks)
                         ipeak = ipeaks(ip);
                         switch ipeak
                             case 1
                                 peak = 'N1';
                             case 2
                                 peak = 'P2';
                         end
                         disp(['Calculating median peak time. Subj = ' num2str(s) ' peak = ' peak ', con = ' num2str(con) ])
                         %calc wave
                         if isstring(electrode)
                            if strcmp(electrode,'GFP')
                                %calculate GFP
                                %wave = ...
                            end
                        elseif numel(electrode)==1
                            wave = nanmean(grandMatrix(:,s,con,:,electrode),4);
                        elseif numel(electrode)>1
                            wave = nanmean(mean(grandMatrix(:,s,con,:,electrode),5),4);
                        end
                         if LP 
                            wave = LPF(wave,srate,20);
                         end
                         subplot(failed_sph{s,con,ipeak})
                         title([peak ' as median'])
                         grandcon_peak_smpls(s,con,ipeak) = round(nanmedian(grandcon_peak_smpls(s,:,ipeak)));
                         grandcon_peak_times(s,con,ipeak) = t(grandcon_peak_smpls(s,con,ipeak));
                         grandcon_peak_amps(s,con,ipeak) = wave(grandcon_peak_smpls(s,con,ipeak));
                         grandcon_as_median(s,con,ipeak) = 1;
                         %plot
                         hold on
                         plot(grandcon_peak_times(s,con,ipeak),grandcon_peak_amps(s,con,ipeak),'o','MarkerSize',10,'Color',Colors{ipeak});
                         text(grandcon_peak_times(s,con,ipeak),grandcon_peak_amps(s,con,ipeak)+0.5,peak,'Color',Colors{ipeak});
                         %now calc and plot prevcons
                         %calc relevant indexes:
                         con_spi = conord(con);%find relevant subplot index:
                         sbplots = size(grandMatrix,3)*(1:(size(grandMatrix,4)-1))+con_spi;
                         sbplots = sbplots(randperm(length(sbplots)));
                         ispp=0;
                         for prevcon = 1:size(grandMatrix,4)
                            %calculate 'wave':
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
                            %calc peak and plot:
                            if ~isnan(wave)
                                if ~(all(~wave))
                                    ispp = ispp+1;
                                    spph = subplot(5,5,sbplots(ispp));
                                    disp(['Elect = ' electrodeName ', subj = ' num2str(s) ', con=' num2str(con) ', prevcon=' num2str(prevcon)])
                                    %calc N1-P2 amp
                                    for ipeak = 1:length(positive_peaks)%run on N1 and P2
                                        smpl = grandcon_peak_smpls(s,con,ipeak);
                                        peak_amp = mean(wave((smpl-ds):(smpl+ds)));                        
                                        %plot wave
                                        hold on;
                                        plot(t,wave);
                                        y = ylim;
                                        rectangle('Position',[t(smpl-ds),y(1), t(smpl+ds)-t(smpl-ds), y(2)-y(1)],'FaceColor',[0.9 0.9 0.9])
                                        plot(t,wave);
                                        xlim([-100 400])
                                        plot([0 0], [y(1) y(2)], '-.k');
                                        plot([-100 400],[0 0],'-.k')
                                        plot([t(smpl-ds) t(smpl+ds)],[peak_amp peak_amp],'Color',Colors{ipeak},'linewidth',2)
                                        peak_smpls(s,con,prevcon,ipeak) = smpl;
                                        if ~isnan(smpl)
                                            peak_amps(s,con,prevcon,ipeak) = peak_amp;
                                            peak_times(s,con,prevcon,ipeak) = t(smpl);
                                        else
                                            peak_amps(s,con,prevcon,ipeak) = nan;
                                            peak_times(s,con,prevcon,ipeak) = nan;
                                        end                            
                                    end

                                    if sbplots(ispp)==25, xlabel('ms'); end

                                end
                            end
                        end
                     end
                 end
            end
            suptitle(['Subject ' num2str(s) '. Block ' block])
            if pauseflag
                disp('Press any key to continue to the next participant')
                pause
            end 
            if ~plotflag
                close
            end
        end
    end

end

