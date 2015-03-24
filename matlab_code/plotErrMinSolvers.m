
clearvars, close all

addpath('..\..\antennaResection\matlab\sharedRoutines.')
%% Constants

experiments=[1 2 3 4 5 6 7] %Index vector on which experiemts should be run
%1 - 4 mics, 4 sounds, 3D,                  calibrated
%2 - 4 mics, 5 sounds, 3D,                  UNcalibrated
%3 - 6 mics, 2 sounds, 2D Mics - 3D Sounds, calibrated
%4 - 8 mics, 2 sounds, 2D Mics - 3D Sounds, UNcalibrated
%5 - 4 mics, 3 sounds, 2D Mics - 3D Sounds, calibrated
%6 - 6 mics, 3 sounds, 3D Mics - 2D Sounds, calibrated
%7 - 8 mics, 3 sounds, 3D Mics - 2D Sounds, UNcalibrated
rigLength=1;
factor=5;                      %Factor that receivers and sender ground truth are scaled up from.

nbrRuns=100; %numer of constellations (i.e. testruns) for each row in nbrMics
stdErr=[10^-7 10^-6 10^-5 10^-4 10^-3];%0.001; %WHite gaussion noise added to theTOA measurements
colors=['b' ;'g' ; 'r'; 'y' ; 'k' ; 'c'; 'm'];

%numer of receivers and transmittors each run


%%
errs=zeros(length(experiments),nbrRuns,length(stdErr));
okChol=zeros(length(experiments),nbrRuns, length(stdErr));

totTimeExp=zeros(length(experiments),1);
hold on

for ii=1:length(experiments)
    
    
    for jj=1:nbrRuns
        
        for kk=1:length(stdErr)
            %CREATING GROUND TRUTH
            %r_gt is receivers and s_gt are senders. Each column is a
            %receiver/sender. Rows are dimensions
            [r_gt, s_gt, D]=makeGroundTruth(experiments(ii),rigLength,factor,stdErr(kk));
            
            %% Run tests
            
            switch ii
                case 1
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 2
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 3
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 4
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 5
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 6
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                case 7
                    tic
                    [r,s,okChol(ii,jj,kk)]=dummyFcn(D);
                    t1=toc;
                otherwise
                    error('Dont ya be hasslin me!')
            end
            
            
            totTimeExp(ii)=totTimeExp(ii)+t1; %NOT over only corect solutions.
            
            %keyboard
            %% solution and gt are only matched up till rotation
            if okChol(ii,jj,kk)
                [~, T]=rigidTransform([r s],[r_gt s_gt]);
                %errsDiffdim(ii,jj)=sqrt(9)*errsDiffdim(ii,jj)/norm([gt.centers gt.coords]); %9* to make it the norm
                
                %TEST so it really is the relative norm
                solAligned=T*[r s; ones(1,size([r s],2))];
                solAligned=solAligned(1:end-1,:);
                errs(ii,jj,kk)=norm(solAligned-[r_gt s_gt],'fro')/norm([gt.centers gt.coords],'fro');
            end
            
            
            % calculate average errors and failure rate
        end
    end
    

    
end

  %to make RMSE
meanRelErr=mean(errs,2);
meanRelErr=reshape(meanRelErr,length(stdErr),length(experiments))
flre=zeros(length(stdErr),length(experiments));
for ii=1:length(stdErr)
    for kk=1:size(nbrRS,1)
        flre(ii,kk)=sum(fail(ii,:,kk))/length(fail(ii,:,kk));
    end
end
fontsize = 14;
set(gca,'fontsize',fontsize);


for ii=1:size(meanRelErr,2)
    loglog(stdErr,meanRelErr(:,ii)',strcat(colors(ii),pStyle(ii,:)),'LineWidth',2);
    hold on;
end
h1=gca;
set(h1,'YAxisLocation','right')
set(h1,'XLim',[stdErr(1)*0.5 stdErr(end)*2])


xlabel('standard deviation or error','fontsize',19);
ylabel('LINES - mean rel. error ','fontsize',19);
%title('Farfield Approximation performance');
legend('6 r - 3 s','7 r - 4 s','20 r - 10 s','Location','NorthWest')

h2 = axes('Position',get(h1,'Position'));
bar(log(stdErr),flre(),0.4);
set(h2,'YAxisLocation','left','Color','none','XTickLabel',[]);
set(h2,'XLim',log(get(h1,'XLim')),'Layer','top');
set(h2,'YLim',[0 0.5],'Layer','top');
set(gca,'fontsize',fontsize);
ylabel('BARS - Complex solution rate','fontsize',19); 
    
    
 %7 - 8 mics, 3 sounds, 3D Mics - 2D Sounds, UNcalibrated  
legend('4 rec - 4 sender, 3D, calib', ...
    '4 rec - 5 sender, 3D,  uncalib', ...
    '6 2D rec - 2 3D sender, calib', ...
    '8 2D rec - 2 3D sender, uncalib', ...
    '4 2D rec - 3 3D sender, calib',...
    '6 3D rec - 3 2D sender, calib',...
    '8 3D rec - 3 2D sender, uncalib')
hold off