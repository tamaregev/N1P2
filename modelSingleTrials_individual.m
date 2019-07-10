function [ SFs, Errs, data, model, sigmas, taus ] = modelSingleTrials_individual( ExpN, Folder, RAdate )
%Compare calculated RA to single trials data

switch ExpN
    case 1
    case 2
    case 3
        Definitions_N1P2
end

whichpeaks = {'N1','P2'};

%load data single trials
Expdir = [Folder filesep 'Exp' num2str(ExpN)];
load([Expdir filesep 'singleAmps'])

%load model - 
loadFolder = [modelFolder RAdate filesep];
fprintf(['Loading...']);tic
load([loadFolder 'RAcat'])
fprintf(['Done in %4.1f sec \n'], toc)
%takes about 1.5 minutes to load RAcat1
load([loadFolder 'Params'])
load([loadFolder 'Metadata'])%artIndss, seqIndss, smplss, stimCodess
load([loadFolder 'metadatacat']);

%reduce size of model :
itaumax = 25;
RAcat = RAcat(:,:,:,:,1:itaumax);
taus = taus(1:itaumax);

%find best scaling factor flags:
bestSFflag = true;
biasFlag = true;%add an intercept to the scaling factor (two scaling params)
if biasFlag
    SFs = nan(whichSubjects(end),length(sigmas),length(taus),length(whichpeaks),2);
else    
    SFs = nan(whichSubjects(end),length(sigmas),length(taus),length(whichpeaks));
end
Errs = nan(whichSubjects(end),length(sigmas),length(taus),length(whichpeaks));

for ipeak = 1:length(whichpeaks)
    disp([whichpeaks(ipeak)])
    for s = whichSubjects
        disp(['Subject ' num2str(s)])
        data = singleAmps(s,:,:,ipeak); % data is subjs x types x timepoints
        data = data(:);
        isArt = isArtefact(s,:,:);isArt = isArt(:);   
        seqi = seqIndcat(s,:,:);seqi=seqi(:);

        for isig = 1:length(sigmas)
           % disp(['sigma = ' num2str(sigmas(isig))])
            for itau = 1:length(taus)
                %RAcat is %subjs x types x timepoints
                model = 1-RAcat(s,:,:,isig,itau);
                model = model(:);

                %scatter plot before cleaning
                if 0
                    x=model; y=data;
                    ERPfigure;
                    plot(x(isArt==1),y(isArt==1),'.','markersize',12,'Color','r')
                    hold on
                    plot(x(seqi<10),y(seqi<10),'.','markersize',16,'Color','y')
                    plot(x(seqi<5),y(seqi<5),'.','markersize',12,'Color','g')
                    plot(x(isArt==0 & seqi>5),y(isArt==0 & seqi>5),'.','markersize',6,'Color','b')
                    title([whichpeaks{ipeak} '. Sigma = ' num2str(sigmas(isig)) ', tau = ' num2str(taus(itau))])
                    xlabel('model');ylabel('data')
               end            
                %exclude artifacts and initial conditions (f for fix):
                modelf = model(isArt==0 & seqi>10); dataf = data(isArt==0 & seqi>10);

                %scatter plot after cleaning
                if 0
                    ERPfigure;
                    x=modelf; y=dataf;
                    plot(x,y,'.')
                    title([whichpeaks{ipeak} '. Sigma = ' num2str(sigmas(isig)) ', tau = ' num2str(taus(itau))])
                    xlabel('model');ylabel('data')
                end

                %model = modelf; data = dataf;
                if bestSFflag
                    if biasFlag
                        newmodel = [ones(size(modelf)),modelf];
                        SF = newmodel\dataf;%performs SS linear regression
                        SFs(s,isig,itau,ipeak,:) = SF;
                        newmodel = [ones(size(modelf)), modelf];
                        predict = newmodel*SF;
                        Errs(s,isig,itau,ipeak) = rms((dataf - predict));
                    else
                        SF = model\dataf;%performs SS linear regression
                        SFs(s,isig,itau,ipeak) = SF;
                        predict = modelf*SF;                    
                        Errs(s,isig,itau,ipeak) = rms((dataf - predict)/mean(ste));    
                    end
                else
    %                 SF = max(abs(data))*sign(data(1));
    %                 SFs(isig,itau,ipeak) = SF;
    %                 %SFs(isig,itau,ipeak) = max(abs(data))*sign(data(1))/max(model);
                end
            end
        end
    end
end

end

