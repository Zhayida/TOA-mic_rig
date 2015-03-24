
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
      
rigLength=1;
factor=5;                      %Factor that receivers and sender ground truth are scaled up from.

nbrRuns=100; %numer of constellations (i.e. testruns) for each row in nbrMics
stdErr=0;%0.001; %WHite gaussion noise added to the TDOA measurements
colors=['b' ;'g' ; 'r'; 'k'];
pointStyle=['d'; 's'; '<'; 'x'];
labels={'2 rigs - 4 send, 3D, calib'; 
        '2 rigs - 5 send, 3D, uncalib',
        '3 rigs in 2D - 2 send in 3D, calib',
        '4 rigs in 2D - 2 send in 3D, uncalib'};

 %numer of receivers and transmittors each run


%%
  errs=zeros(max(experiments),nbrRuns);
okChol=zeros(max(experiments),nbrRuns);

totTimeExp=zeros(max(experiments),1);
hold on
 
for ii=experiments
    
    
    for jj=1:nbrRuns
        
        currExp=ii %FOr display info
        currRun=jj %FOr display info
        %CREATING GROUND TRUTH
        %r_gt is receivers and s_gt are senders. Each column is a
        %receiver/sender. Rows are dimensions
        [r_gt, s_gt, D]=makeGroundTruth(ii,rigLength,factor,0);
        
        %% Run tests
        
        switch ii
            case 1
                
                data.d=D;
                data.l=rigLength;
                
                tic
                [rpos,spos,okChol(ii,jj)] = solver_4_4(data);
                t1=toc;
                
                bestErr=Inf;
                for kk=1:size(rpos,2)
                    [~, ~,solTransformed]=rigidTransform([rpos{kk} spos{kk}],[r_gt s_gt]);
                    currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                    
                    if currErr < bestErr
                        bestErr=currErr;
                        errs(ii,jj)=currErr;
                    end
                end
                
            case 2
                data.d=D;
                data.l=[]; %Just a precaution. kill old value                
                tic
                [rpos,spos,okChol(ii,jj)] = solver_4_5(data);
                t1=toc;
                
                bestErr=Inf;
                for kk=1:size(rpos,2)
                    [~, ~,solTransformed]=rigidTransform([rpos{kk} spos{kk}],[r_gt s_gt]);
                    currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                    
                    if currErr < bestErr
                        bestErr=currErr;
                        errs(ii,jj)=currErr;
                    end
                end
                
            case 3
                tic
                [sols]=toa_rig_6_2(D,[rigLength rigLength rigLength]);
                t1=toc;
                
                okChol(ii,jj)=0;
                bestErr=Inf;
                for kk=1:size(sols.R,2)
                    
                    if isreal([sols.R{kk} sols.S{kk}])
                        okChol(ii,jj)=1;
                        [~, ~,solTransformed]=rigidTransform([sols.R{kk} sols.S{kk}],[r_gt s_gt]);
                        currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                        
                        if currErr < bestErr
                            bestErr=currErr;
                            errs(ii,jj)=currErr;
                        end
                    end
                end
                
                if bestErr >10^-5
                   keyboard 
                end
            case 4
                tic
                clear sols
                [sols]=toa_rig_8_2(D);
                t1=toc;
                
                okChol(ii,jj)=0;
                bestErr=Inf;
                for kk=1:size(sols.R,2)
                    
                    if isreal([sols.R{kk} sols.S{kk}])
                        okChol(ii,jj)=1;
                        [~, ~,solTransformed]=rigidTransform([sols.R{kk} sols.S{kk}],[r_gt s_gt]);
                        currErr=norm(solTransformed-[r_gt s_gt],'fro')/norm([r_gt s_gt],'fro');
                        
                        if currErr < bestErr
                            bestErr=currErr;
                            errs(ii,jj)=currErr;
                        end
                    end
                end
            case 5
                tic
                [r,s,okChol(ii,jj)]=dummyFcn(D);
                t1=toc;
            case 6
                tic
                [r,s,okChol(ii,jj)]=dummyFcn(D);
                t1=toc;
            case 7
                tic
                [r,s,okChol(ii,jj)]=dummyFcn(D);
                t1=toc;
            otherwise
                error('Dont ya be hasslin me!')
        end
        
        
        totTimeExp(ii)=totTimeExp(ii)+t1; %NOT over only corect solutions.
        

        
    end
    
    
    
    
end

% hold on
% for ii=experiments
%     nbin = 10;
%     linewidth = 3;
%     fontsize = 14;
%     toPlot=errs(ii,errs(ii,:)~=0);
%     [aa,bb] = hist(log10(errs(ii,:)),nbin);
%     plot(bb,aa,strcat(colors(ii),'-'),'linewidth',linewidth);
%     set(gca,'fontsize',fontsize);
%     xlabel('log_{10}(rel errors)','fontsize',19);
%     ylabel('Frequency','fontsize',19);
% end
% hold off

maxbin=max(log10(errs(errs~=0)));
minbin=min(log10(errs(errs~=0)));
figure; hold on
for ii=experiments
     nbins=20
    linewidth = 2;
    fontsize = 14;
    toPlot=errs(ii,errs(ii,:)~=0);
    [aa,bb] = hist(log10(toPlot),linspace(minbin,maxbin,nbins));
    plot(bb,aa,strcat(colors(ii),'-',pointStyle(ii)),'linewidth',linewidth);
    set(gca,'fontsize',fontsize);
    %set(gca,'XTick',[-16:1:-6])
    xlabel('log_{10}(rel errors)','fontsize',19);
    ylabel('Frequency','fontsize',19); 
end
axis([-16 0,0,50])
legend(labels{experiments})
hold off