function EV = ExpVar(x,y,SF,ws)
%x=modelf; y=dataf;
%plot(x,y,'.')
%ws=1025;
[sm,I]=sort(x);
sr=y(I);
b=dpss(ws,2,1);b=b/sum(b);
fsr=filtfilt(b,1,sr);
%plot(sm,fsr,'.')
ssm=sm*SF(2)+SF(1);
efsr=fsr(ceil(ws/2):end-floor(ws/2));
essm=ssm(ceil(ws/2):end-floor(ws/2));
EV = 1-mean((efsr-essm).^2)/var(efsr);
end