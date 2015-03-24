
clearvars, close all

addpath('..\..\antennaResection\matlab\sharedRoutines')  %Different depending where you put antennaresection and micrig
addpath('.\solvers')
addpath('C:\Users\simonb\Documents\multipol')
%% Constants

experiments=[1 2 3 4]; %[1 2 3];% 2 3 4];% 5 6 7] %Index vector on which experiemts should be run
                                %1 - 4 mics, 4 sounds, 3D,*                  calibrated    WORKS         
                                %2 - 4 mics, 5 sounds, 3D,*                  UNcalibrated WORKS most of cases. SOmetimes (1/30) it only finds complex solutions.
                                %3 - 6 mics, 2 sounds, 2D Mics - 3D Sounds,* calibrated   WORKS    
                                %4 - 8 mics, 2 sounds, 2D Mics - 3D Sounds,* UNcalibrated WORKSd 
                                %5 - 4 mics, 3 sounds, 2D Mics - 3D Sounds,* calibrated    
                                %6 - 6 mics, 3 sounds, 3D Mics - 2D Sounds*, calibrated 
                                %7 - 8 mics, 3 sounds, 3D Mics - 2D Sounds*, UNcalibrated
      
rigLength=0.2;
factor=1;                      %Factor that receivers and sender ground truth are scaled up from.
stdDevNoiseVec=10.^[-7 -6 -5 -4]
nbrRuns=100; %numer of constellations (i.e. testruns) for each row in nbrMics
%stdErr=0;%0.001; %WHite gaussion noise added to the TDOA measurements
colors=['b' ;'g' ; 'r'; 'k'];
pStyle=['d-'; 's-'; '<-'; 'x-'];
labels={'2 rigs - 4 send, 3D, calib'; 
        '2 rigs - 5 send, 3D, uncalib',
        '3 rigs in 2D - 2 send in 3D, calib',
        '4 rigs in 2D - 2 send in 3D, uncalib'};

 %numer of receivers and transmittors each run


%%
  errs=zeros(length(stdDevNoiseVec),max(experiments),nbrRuns);
okChol=zeros(length(stdDevNoiseVec),max(experiments),nbrRuns);

totTimeExp=zeros(max(experiments),1);
hold on
for hh=1:length(stdDevNoiseVec)
    for ii=experiments
        
        
        for jj=1:nbrRuns
            
            currExp=ii %FOr display info
            currRun=jj %FOr display info
            %CREATING GROUND TRUTH
            %r_gt is receivers and s_gt are senders. Each column is a
            %receiver/sender. Rows are dimensions
            [r_gt, s_gt, D]=makeGroundTruth(ii,rigLength,factor,stdDevNoiseVec(hh));
            
            %% Run tests
            
            switch ii
                case 1
                    
                    data.d=D;
                    data.l=rigLength;
                    
                    tic
                    [rpos,spos,okChol(hh,ii,jj)] = solver_4_4(data);
                    t1=toc;
                    
                    bestErr=Inf;
                    for kk=1:size(rpos,2)
                        [~, ~,solTransformed]=rigidTransform([rpos{kk} spos{kk}],[r_gt s_gt]);
                        currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                        
                        if currErr < bestErr
                            bestErr=currErr;
                            errs(hh,ii,jj)=currErr;
                        end
                    end
                    
                case 2
                    data.d=D;
                    data.l=[]; %Just a precaution. kill old value
                    tic
                    [rpos,spos,okChol(hh,ii,jj)] = solver_4_5(data);
                    t1=toc;
                    
                    bestErr=Inf;
                    for kk=1:size(rpos,2)
                        [~, ~,solTransformed]=rigidTransform([rpos{kk} spos{kk}],[r_gt s_gt]);
                        currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                        
                        if currErr < bestErr
                            bestErr=currErr;
                            errs(hh,ii,jj)=currErr;
                        end
                    end
                    
                case 3
                    tic
                    [sols]=toa_rig_6_2(D,[rigLength rigLength rigLength]);
                    t1=toc;
                    
                    okChol(hh,ii,jj)=0;
                    bestErr=Inf;
                    for kk=1:size(sols.R,2)
                        
                        if isreal([sols.R{kk} sols.S{kk}])
                            okChol(hh,ii,jj)=1;
                            [~, ~,solTransformed]=rigidTransform([sols.R{kk} sols.S{kk}],[r_gt s_gt]);
                            currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                            
                            if currErr < bestErr
                                bestErr=currErr;
                                errs(hh,ii,jj)=currErr;
                            end
                        end
                    end
                    
 
                case 4
                    tic
                    clear sols
                    [sols]=toa_rig_8_2(D);
                    t1=toc;
                    
                    okChol(hh,ii,jj)=0;
                    bestErr=Inf;
                    for kk=1:size(sols.R,2)
                        
                        if isreal([sols.R{kk} sols.S{kk}])
                            okChol(hh,ii,jj)=1;
                            [~, ~,solTransformed]=rigidTransform([sols.R{kk} sols.S{kk}],[r_gt s_gt]);
                            currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                            
                            if currErr < bestErr
                                bestErr=currErr;
                                errs(hh,ii,jj)=currErr;
                            end
                        end
                    end
                case 5
                    tic
                    [r,s,okChol(hh,ii,jj)]=dummyFcn(D);
                    t1=toc;
                case 6
                    tic
                    [r,s,okChol(hh,ii,jj)]=dummyFcn(D);
                    t1=toc;
                case 7
                    tic
                    [r,s,okChol(hh,ii,jj)]=dummyFcn(D);
                    t1=toc;
                otherwise
                    error('Dont ya be hasslin me!')
            end
            
            
            totTimeExp(ii)=totTimeExp(ii)+t1; %NOT over only corect solutions.
            
            
            
        end
        
        
        
        
    end
    
end

%% Plotting
figure
log10(errs);
%errs=zeros(length(stdDevNoiseVec),max(experiments),nbrRuns);
meanlog10Errs=zeros(length(stdDevNoiseVec) ,max(experiments));
for ii=1:length(stdDevNoiseVec) 
    for jj=1:max(experiments)
        errsValid=errs(ii,jj,:);
        errsValid=errsValid(errsValid~=0); %take out the ones ~=0.
        meanlog10Errs(ii,jj)=mean(log10(errsValid));
    end
end

clear ii jj  hh; 
flre=zeros(length(stdDevNoiseVec) ,max(experiments));
for ii=1:length(stdDevNoiseVec)
    for jj=1:max(experiments)
        flre(ii,jj)=sum(errs(ii,jj,:)==0)/nbrRuns;
    end
end


fontsize = 14;
set(gca,'fontsize',fontsize);
meanErrs=(meanlog10Errs); %Just to get the ewaxis right.
for ii=1:size(meanlog10Errs,2)
    plot(log10(stdDevNoiseVec),meanErrs(:,ii),strcat(colors(ii),pStyle(ii,:)),'LineWidth',2);
    hold on;
end



h1=gca;
set(h1,'YAxisLocation','right')
set(h1,'XLim',log10([stdDevNoiseVec(1)*0.3 stdDevNoiseVec(end)*3]))
set(h1,'YLim',[-6 -1])
set(gca,'XTick',[-7:1:-4])
set(gca,'YTick',[-6:1:-1])

xlabel('log10(standard deviation)','fontsize',19);
ylabel('LINES - log10(rel. error) ','fontsize',19);
%title('Farfield Approximation performance');
%legend('6 r - 3 s','7 r - 4 s','20 r - 10 s','Location','NorthWest')

h2 = axes('Position',get(h1,'Position'));
colormap([0 0 1; 0 1 0 ; 1 0 0  ; 0 0 0])
bar(log10(stdDevNoiseVec),flre(),0.9);
set(h2,'YAxisLocation','left','Color','none','XTickLabel',[]);
set(h2,'XLim',get(h1,'XLim'),'Layer','top');
set(h2,'YLim',[0 0.3],'Layer','top');
set(gca,'fontsize',fontsize);
set(gca,'YTick',[0:0.1:0.3])
ylabel('BARS - Complex solution rate','fontsize',19); 

