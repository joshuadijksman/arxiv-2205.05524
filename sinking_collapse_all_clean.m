




% some details regarding conversions etc 
D2si = 0.01^2/60; % convert creep slopes fitted to SI units (measurements are in cm/min)
Eta2si = 0.01/60; % convert linear slopes fitted to SI units (measurements are in cm/min)

w2pa_30 = 9.81*0.001/(3.1415*0.015^2); % surface area of plunger 30mm cyl
w2pa_35 = 9.81*0.001/(3.1415*0.0175^2); % surface area of plunger 35mm cyl
w2pa_23 = 9.81*0.001/(3.1415*0.0115^2); % surface area of sphere 23mm



%% read data - keep script in same directory as the data files

fnames = dir('*.dat');
numfiles = length(fnames);

densarr = NaN(numfiles,1);
sizearr = NaN(numfiles,1);
addedweightarr = NaN(numfiles,1);
lengarr = NaN(numfiles,1);

timearr = NaN(numfiles, 50); % make enough space
displarr = NaN(numfiles, 50);
errarr = NaN(numfiles, 50);

% for creep curves
slopearr = NaN(numfiles,1);
offsetarr = NaN(numfiles,1);
errorarr = NaN(numfiles,2);

% for linear curves
slopearr_lin = NaN(numfiles,1);
offsetarr_lin = NaN(numfiles,1);
errorarr_lin = NaN(numfiles,2);
conf_lin = NaN(numfiles,2,2);

%settings to extract best fit exponent
bexpy = 0.2; % lower exp
expyr = 0.025; % step exponent
eexpy = 1.4; % max exp
expyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % save all exp fit rsquares here
stdarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % stats of residuals
kurtarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % kurtosis stats of residuals
skewarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % skewness stats of residuals
distopyarr = NaN(numfiles,1+round((eexpy-bexpy)/expyr)); % distance to optimum

%save exponents to memory
bestexparr = NaN(numfiles,1); 
optexparr = NaN(numfiles,1);


typearr = []; % save type: cyl or sphere


%run over all data files that contain info of individual creep experiments
for i=1:numfiles
    
    if ~mod(i,10), fprintf('numfile # %d\n', i); end
    
    % units all in centimeters and minutes!!!!!!!!!!
    readfile = char(fnames(i).name);
    
    densarr(i) = str2num(readfile(1:4))/1000;
    typearr = [typearr, readfile(5)];
    sizearr(i) = str2num(readfile(6:7));
    addedweightarr(i) = str2num(readfile(8:9));
    
    tmp = importdata(readfile);
    tmpdata = tmp.data;
    tmpstring = tmp.textdata;
    tmpstring = tmpstring{1};
    datl = size(tmpdata);
    
    lengarr(i) = datl(1);
    timearr(i,1:datl) = tmpdata(:,1);
    displarr(i,1:datl) = tmpdata(:,2);
    

    if contains(tmpstring, 'err') % if there is an error estimate, it will be in the third column.
        errarr(i,1:datl(1)) = tmpdata(:,3);
    end
    
    %%%%%%%%%%%% we first assume square root behavior:
    % delta = sqrt(Dx)
    % delta^2 = Dx
    % x = time
    % delta = displacement

    xtmp = [timearr(i,1:lengarr(i))];
    deltatmp = [displarr(i,1:lengarr(i)).^2]; % square first then do linear fit. Prefactor is D
    
    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');
    
    tmp = coeffvalues(fitobj);    
    %prefactor is slope on y^2(t) plot
    slopearr(i) = tmp(1);
    offsetarr(i) = tmp(2);
    errorarr(i,1) = gof2.rmse;%/datl(1);
    errorarr(i,2) = gof2.rsquare;
    
    %%%%%%%%%%%%%%%% now fit the linear parts
    deltatmp = [displarr(i,1:lengarr(i))];
    
    % fit poly1, offset should be zero hence added point at origin
    [fitobj,gof2] = fit(xtmp',deltatmp','poly1');
    
    %figure
    %plot(fitobj,xtmp',deltatmp')
    
    tmp = coeffvalues(fitobj);    
    slopearr_lin(i) = tmp(1);
    offsetarr_lin(i) = tmp(2);
    errorarr_lin(i,1) = gof2.rmse;%/datl(1);
    errorarr_lin(i,2) = gof2.rsquare;
    conf_lin(i,:,:) = confint(fitobj,0.95); 
    
    %%%%%%%%%%%%%%%% fit a range of different exponents
    xtmp = timearr(i,2:lengarr(i));
    deltatmp = displarr(i,2:lengarr(i));
    ctr = 0;
    for expy = bexpy:expyr:eexpy
        ctr = ctr + 1;
        opty = fitoptions('power2','Lower', [0,expy,-10], 'Upper', [Inf,expy,10]);
       [fitobj,gof2] = fit(xtmp',deltatmp','power2',opty);
       expyarr(i,ctr) = gof2.rsquare;
       %thought: check stats of residuals
       resi = fitobj(xtmp)-deltatmp';
       stdarr(i,ctr) = nanstd(resi);
       kurtarr(i,ctr) = kurtosis(resi);
       skewarr(i,ctr) = skewness(resi);  
    end
    [~,indexy] = max(expyarr(i,:));
    bestexparr(i) = bexpy+(indexy-1)*expyr;
    
    distopyarr(i,:) = sqrt(skewarr(i,:).^2 + (kurtarr(i,:)-3).^2 + (gradient(expyarr(i,:)).^2));
    [~,indexy] = min(distopyarr(i,:));
    optexparr(i) = bexpy+(indexy-1)*expyr;
end

%% extract locations of cylinder and sphere

maskcyl30 = find(sizearr == 30);
maskcyl35 = find(sizearr == 35);
masksp = find(sizearr == 23);
masksv = find(sizearr == 24); % new sieved particles were done with 23mm sphere though. 24mm is just for labeling

densarruniq = unique(densarr); % set of all densities probed. Be careful: 7.02 exists for sphere and cylinder!

addweiuniq = unique(addedweightarr);

mindens = min(unique(densarr));
maxdens = max(unique(densarr));
densrange = maxdens-mindens;

prefacarr = NaN(length(unique(densarr)),1); % to be filled manually with rescaling factors


%% calculate stresses per experiment

denscyl = 1190; %acrylic cylinder
denssphere = 913;  %polypropylene sphere
denswater = 1006;  % at standard T, p. Ignoring that we use Oxford tap water.
weightrodsp = 0.00224; % rod for spheres is longer than for cylinders
weightrodcyl = 0.00176; % also in kg
weighttray = 0.0134;  % also in kg
weightplunger30 = (3.1415*0.015^2)*0.02*(denscyl-denswater); %compute effective weight from delta dens: pi*r^2*h*Delta-rho
weightplunger35 = (3.1415*0.0175^2)*0.02*(denscyl-denswater);
weightplunger23 = (4*3.1415*0.0115^3)*0.333*(denssphere-denswater); %(4/3)*pi*r^3*Delta-rho
surfplung30 = 3.1415*0.015^2; % surface area of plunger 30mm cyl
surfplung35 = 3.1415*0.0175^2; % surface area of plunger 35mm cyl
surfplung23 = 3.1415*0.0115^2; % surface area of sphere 23mm

stressplung30 = 9.81*(weightplunger30+weightrodcyl+weighttray)/surfplung30;
stressplung35 = 9.81*(weightplunger35+weightrodcyl+weighttray)/surfplung35;
stressplung23 = 9.81*(weightplunger23+weightrodsp+weighttray)/surfplung23;
% apparently the largest densities have the cylinder only partially
% immersed so the stress calculation should be different.

stressarr(maskcyl30) = 9.81*(weightplunger30+weightrodcyl+weighttray + addedweightarr(maskcyl30)/1000)/surfplung30;
stressarr(maskcyl35) = 9.81*(weightplunger35+weightrodcyl+weighttray + addedweightarr(maskcyl35)/1000)/surfplung35;
stressarr(masksp) = 9.81*(weightplunger23+weightrodsp+weighttray + addedweightarr(masksp)/1000)/surfplung23;
stressarr(masksv) = 9.81*(weightplunger23+weightrodsp+weighttray + addedweightarr(masksv)/1000)/surfplung23;


% data interpretation settings for all figures - indicate the edges of the
% two regimes. These numbers are determined manually
crtreshsp = 6.82;
lintreshsp = 7.15;

platysp = 1.6E-5;
sloperhosp = 13;
sloperhocyl = 26;
platycyl = 2E-6;
Dsigmsp = 3.75E-5;
Dsigmcyl = 9E-4;

%% make the first overview figure 


figure(1)

subplot 222
cmap = colormap(parula(length(unique(densarr(masksp)))));

for i=1:4:length(masksp)
    
    indexer = find(unique(densarr(masksp))==densarr(masksp(i)));
    
    if densarr(masksp(i)) < lintreshsp
       
        scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k'); %,'MarkerEdgeColor',cmap(indexer,:))
        hold on
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k')%,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksp(i))) 
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'k')%,'Color',cmap(indexer,:))
            fprintf('errorrrrrrrrrr %d\n', expy)
        end
        
        
    end
    
end


box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(3250,.015,'(b)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([0 4000])%xlim([20,6000])
ylim([0,0.1])%ylim([3E-3,.15])
ax=gca;
ax.FontSize = 14;
hold off



subplot 224

for i=length(masksp):-1:1
    
    indexer = find(unique(densarr(masksp))==densarr(masksp(i)));
    
    if densarr(masksp(i)) < lintreshsp
       
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-ok')

        %convert D to proper units, but also the war time and displacement data
        hold on
        
    end
    
end
plot([10,4000],[10,4000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
text(3250,500,'(c)','FontSize',18)
ax=gca;
ax.FontSize = 14;
hold off

% save
%set(gcf, 'Color', 'None') %overlapping figure in APS compiler were not
%solved with this line
print('fig1', '-depsc','-r600')
print('fig1', '-dpdf', '-bestfit')

%%

figure(2)

% first plot the D dynamics for both cases

subplot(2,2,[1 2])


plot([6.3,crtreshsp],[platysp,platysp], '-.k') % for spheres
hold on
plot([crtreshsp,crtreshsp],[1E-3,1E3], ':k')


for fixmass  = 0:0 % 1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            psh1 = scatter([densarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],'b','Fill')
            %convert D to proper units, but also the war time and displacement data
            hold on
        end

    end
    
    xtmp = crtreshsp:0.01:7.2;
    plot(xtmp,platysp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')

end

% also for the cylinder
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')


plot([6.3,crtreshsp],[platycyl,platycyl], '-.k') % for cyl

for fixmass  = 0:0
 
    for i=1:length(maskcyl30)

        if addedweightarr(maskcyl30(i)) == fixmass,
            pcyl1 = scatter([densarr(maskcyl30(i))],[slopearr(maskcyl30(i))*D2si],[30],'rs', 'Fill')

            %convert D to proper units, but also the war time and displacement data
            hold on
        end

    end
    
    xtmp = crtreshsp:0.01:8.1;
    plot(xtmp,platycyl*exp(-sloperhocyl*(xtmp-crtreshsp)),'-k')


 

end
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')

xlim([6.51,7.17])
ylim([1e-8,5E-5])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D [m^2/s]')
box on
legend([psh1,pcyl1], 'sphere','cylinder','Location','southwest')


text(7.0,1E-5,'\sigma = 355 Pa','Color','blue','FontSize',14)
text(6.6,7E-7,'\sigma = 247 Pa','Color','red','FontSize',14)


% save
print('fig2', '-depsc','-r600')
print('fig2', '-dpdf', '-bestfit')


%% rescaled plot, first sphere

figure(3)

% first sphere data

bluey = subplot(2,2,1);

xtmp = 300:1:600;
plot(xtmp,0.65E-10*exp(xtmp./35),'-.k')
hold on

for i=1:length(masksp)

    if densarr(masksp(i)) < lintreshsp, %== 6.96,
        scatter([stressarr(masksp(i))],[slopearr(masksp(i))*D2si],[30],[densarr(masksp(i))],'Fill')

        hold on
    end

end

xlim([340,510])
ylim([5e-7,1E-4])
set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
text(485,1E-6,'(a)','FontSize',18)
box on
colorbar('EastOutside')
colorbar('Ticks',[6.8,6.9,7.0],...
         'TickLabels',{'6.8','6.9','7.0'})
%caxis([350,450])
ax=gca;
ax.FontSize = 14;

bluey = subplot(2,2,3);

plot([6.3,crtreshsp],[Dsigmsp,Dsigmsp], '-.k'); 
hold on
plot([crtreshsp,crtreshsp],[1E-9,1E-3], ':k')
xtmp = crtreshsp:0.01:7.2;
plot(xtmp,Dsigmsp*exp(-sloperhosp*(xtmp-crtreshsp)),'-k')


for fixmass  = 0:1:60
 
    for i=1:length(masksp)

        if addedweightarr(masksp(i)) == fixmass,
            scatter([densarr(masksp(i))],[slopearr(masksp(i))*D2si]/(platysp*exp(stressarr(masksp(i))/35)),[30],[stressarr(masksp(i))],'Fill')

            hold on
        end

    end

end
xlim([6.7,7.17])
ylim([1e-7,1E-4])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D/D_0 exp(\sigma/\sigma_s) [-]')
text(7.1,4E-5,'(c)','FontSize',18)
box on
colorbar('EastOutside')
caxis([350,450])
ax=gca;
ax.FontSize = 14;

% CYL 30
 
redy = subplot(2,2,2);
xtmp = 200:1:530;
plot(xtmp,3E-11*exp(xtmp./35),'-.k')
hold on

for i=1:length(maskcyl30)

    if densarr(maskcyl30(i)) < lintreshsp, %== 6.96,
        scatter([stressarr(maskcyl30(i))],[slopearr(maskcyl30(i))*D2si],[30],[densarr(maskcyl30(i))],'s', 'Fill')

        hold on
    end

end

xlim([230,512])
ylim([1e-8,2E-5])
set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('D [m^2/s]')
text(470,2.5E-8,'(b)','FontSize',18)
box on
colorbar('EastOutside')
%caxis([350,450])
ax=gca;
ax.FontSize = 14;

redy = subplot(2,2,4);
plot([6.3,crtreshsp],[Dsigmcyl,Dsigmcyl], '-.k'); 
hold on
xtmp = crtreshsp:0.01:8.1;
plot(xtmp,Dsigmcyl*exp(-sloperhocyl*(xtmp-crtreshsp)),'-k')


for fixmass  = 0:1:60
 
    for i=1:length(maskcyl30)

        if addedweightarr(maskcyl30(i)) == fixmass,
            scatter([densarr(maskcyl30(i))],[slopearr(maskcyl30(i))*D2si]/(platycyl*exp(stressarr(maskcyl30(i))/35)),[30],[stressarr(maskcyl30(i))],'s', 'Fill')

            %convert D to proper units, but also the war time and displacement data
            hold on
        end

    end
    
end

plot([crtreshsp,crtreshsp],[1E-9,1E-1], ':k')

xlim([6.51,7.17])
ylim([5e-7,5E-3])
set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('D/D_0 exp(\sigma/\sigma_s) [-]')
text(7.07,1.5E-3,'(d)','FontSize',18)
box on
colorbar('EastOutside')
caxis([200,450])
ax=gca;
ax.FontSize = 14;

MapRed = [linspace(0,1,256).^1',zeros(256,2)];
MapBlue = [zeros(256,2),linspace(0,1,256)'];
colormap(redy, MapRed);
colormap(bluey, MapBlue);

% save
print('fig3', '-depsc','-r600')
print('fig3', '-dpdf', '-bestfit')




%% linear sinking part

figure(4)

subplot 221

%cmap = colormap(parula(length(densarruniq)));

cmap = colormap(parula(length(unique(densarr(masksp)))));

for i=1:1:length(masksp)
    
    indexer = find(unique(densarr(masksp))==densarr(masksp(i)));

    %if densarr(masksp(i)) > lintreshsp && densarr(masksp(i)) < 7.45
    if addedweightarr(masksp(i)) == 0 && densarr(masksp(i)) > 7.1
        
        
            densarr(masksp(i));
               
        scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,'MarkerEdgeColor',cmap(indexer,:))
        hold on
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksp(i)))
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'Color',cmap(indexer,:))
        end
        hold on
        
    end
    
end
%caxis([lintreshsp,7.4])

box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(3350,0.003,'(a)','FontSize',18)
ax=gca;
ax.FontSize = 14;



% ETA dynamics in inset

subplot 223

%axes('Position',[.171,.62,.14,.27])

for fixmass  = 0:0
 
    for i=1:length(masksp)
        
        indexer = find(unique(densarr(masksp))==densarr(masksp(i)));

        if addedweightarr(masksp(i)) == fixmass,
            
            fgrav = 9.81*(weightplunger23+weightrodsp+weighttray + addedweightarr(masksp(i))/1000);
                        
            scatter([densarr(masksp(i))],[fgrav/(6*3.1415*0.0115*slopearr_lin(masksp(i))*Eta2si)]/1E6,[50],cmap(indexer,:),'Fill')

            hold on
                     
        end

    end
    

end

xtmp = lintreshsp:0.01:7.5;
ytmp = .02./(7.5-xtmp).^1;
plot(xtmp,ytmp,'-.k')

xlim([lintreshsp,7.5])
ylim([0,0.6])
%set(gca,'Yscale','log')
xlabel('\rho [g/L]')
ylabel('\eta_{eff} [MPa.s]')
xticks([7.2,7.3,7.4])
box on
%ax=gca;
%ax.FontSize = 8;
text(7.2,0.5,'(b)','FontSize',18)
ax=gca;
ax.FontSize = 14;


hold off


% v_s (sigma) dynamics
subplot(2,2,[2,4])

for i=1:length(masksp)

    if densarr(masksp(i)) == 7.45,

        scatter([stressarr(masksp(i))],[slopearr_lin(masksp(i))*Eta2si],[70],stressarr(masksp(i)),'Fill')
        
        hold on
        
        poserr = (slopearr_lin(masksp(i)) - conf_lin(masksp(i),1,1))*Eta2si;
        negerr = (conf_lin(masksp(i),2,1) - slopearr_lin(masksp(i)))*Eta2si;

        errorbar([stressarr(masksp(i))],[slopearr_lin(masksp(i))*Eta2si],...
            [poserr],[negerr],'k')
    end

end

xtmp = 200:600;
plot(xtmp,5.5E-11*exp(xtmp/35),'-.k')
    
box on
xlim([350,510])
ylim([1e-6,1.5E-4])
text(485,1.8E-6,'(c)','FontSize',18)
set(gca,'Yscale','log')
xlabel('\sigma [Pa]')
ylabel('v_s [m/s]')
colorbar('EastOutside')
ax=gca;
ax.FontSize = 14;


hold off



% save
print('fig4', '-depsc','-r600')
print('fig4', '-dpdf', '-bestfit')

%% SI figure log-log scale version Fig 1


figure(5)

subplot 221
cmap = colormap(parula(length(unique(densarr(masksp)))));

for i=1:4:length(masksp)
    
    %indexer = find(densarruniq==densarr(masksp(i)));
    indexer = find(unique(densarr(masksp))==densarr(masksp(i)));
    
    if densarr(masksp(i)) < lintreshsp
       
        scatter(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k'); %,'MarkerEdgeColor',cmap(indexer,:))
        hold on
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01, 'k')%,'Color',cmap(indexer,:))

        if ~isnan(errarr(masksp(i))) 
            errorbar(timearr(masksp(i),:)*60,displarr(masksp(i),:)*0.01,errarr(masksp(i),:)*0.01,'k')%,'Color',cmap(indexer,:))
            fprintf('errorrrrrrrrrr %d\n', expy)
        end
        
        
    end
    
end

set(gca,'Xscale','log')
set(gca,'Yscale','log')

plot([10,100000],0.0015*[10,100000].^0.5,'-.k')
plot([10,100000],0.0015*[10,100000].^0.4,':k')
plot([10,100000],0.0015*[10,100000].^0.6,':k')


box on
xlabel('t [sec]')
ylabel('\delta [m]')
text(2000,.0035,'(a)','FontSize',18)
%set(gca,'Yscale','log')
%set(gca,'Xscale','log')
xlim([20, 5000])
ylim([0.002,0.2])
ax=gca;
ax.FontSize = 14;
hold off



subplot 222

for i=length(masksp):-1:1
    
    indexer = find(unique(densarr(masksp))==densarr(masksp(i)));
    
    if densarr(masksp(i)) < lintreshsp
       
        plot(timearr(masksp(i),:)*60,displarr(masksp(i),:).^2/(10000*slopearr(masksp(i))*D2si),'-ok')

        %convert D to proper units, but also the war time and displacement data
        hold on
        
    end
    
end
plot([10,10000],[10,10000],'-.k')
xlabel('t [sec]')
ylabel('\delta^2/D [sec]')
text(2000,8.5,'(b)','FontSize',18)
xlim([20, 5000])
ylim([4,5000])
ax=gca;
ax.FontSize = 14;
hold off

set(gca,'Xscale','log')
set(gca,'Yscale','log')

% save
print('figSI', '-depsc','-r600')
print('figSI', '-dpdf', '-bestfit')

