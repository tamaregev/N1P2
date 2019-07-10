function [ SFs, Errs, Evs, sigmas, taus, ws ] = modelSingleTrials( ExpN, Folder, RAdate )
%Compare calculated RA to single trials data
ws=2048;
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
%takes about 1.5 minutes to load RAcat
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
    SFs = nan(length(sigmas),length(taus),length(whichpeaks),2);
else    
    SFs = nan(length(sigmas),length(taus),length(whichpeaks));
end
Errs = nan(length(sigmas),length(taus),length(whichpeaks));
Evs = nan(length(sigmas),length(taus),length(whichpeaks));

isArt = isArtefact(whichSubjects,:,:);isArt = isArt(:);   
seqi = seqIndcat(whichSubjects,:,:);seqi=seqi(:);

for ipeak = 1:length(whichpeaks)
    disp([whichpeaks(ipeak)])
    for isig = 1:length(sigmas)
        disp(['sigma = ' num2str(sigmas(isig))])
        for itau = 1:length(taus)
            %RAcat is %subjs x types x timepoints
           data = singleAmps(whichSubjects,:,:,ipeak); % data is subjs x types x timepoints
           data = data(:);
 
            model = 1-RAcat(whichSubjects,:,:,isig,itau);
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
            
            model = modelf; data = dataf;
            if bestSFflag
                if biasFlag
                    newmodel = [ones(size(model)),model];
                    SF = newmodel\data;%performs SS linear regression
                    SFs(isig,itau,ipeak,:) = SF;
                    newmodel = [ones(size(model)), model];
                    predict = newmodel*SF;
                    Errs(isig,itau,ipeak) = rms((data - predict));
                    Evs(isig,itau,ipeak) = ExpVar(model,data,SF,ws);
                else
                    SF = model\data;%performs SS linear regression
                    SFs(isig,itau,ipeak) = SF;
                    predict = model*SF;                    
                    Errs(isig,itau,ipeak) = rms((data - predict)/mean(ste));    
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

