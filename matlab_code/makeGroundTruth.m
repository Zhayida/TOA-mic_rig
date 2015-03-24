function [r_gt, s_gt,D]=makeGroundTruth(experimentType,rigLength, factor,stdErr)
%Makes groudn truth for seven different minimal cases of the dual receiver
%rig calibration problem 
%
%ExperimentType can be a number between 1-7. 
                                %1 - 4 mics, 4 sounds, 3D,                  calibrated
                                %2 - 4 mics, 5 sounds, 3D,                  UNcalibrated
                                %3 - 6 mics, 2 sounds, 2D Mics - 3D Sounds, calibrated 
                                %4 - 8 mics, 2 sounds, 2D Mics - 3D Sounds, UNcalibrated 
                                %5 - 4 mics, 3 sounds, 2D Mics - 3D Sounds, calibrated 
                                %6 - 6 mics, 3 sounds, 3D Mics - 2D Sounds, calibrated 
                                %7 - 8 mics, 3 sounds, 3D Mics - 2D Sounds, UNcalibrated 
    switch experimentType
        case 1
            r_gt=factor*(rand(3,4) -0.5); 
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
                z= 2*rand -1;
                r=sqrt(1-z^2);
                rigMove=rigLength*[r*cos(phi) ; r*sin(phi); z];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            s_gt=factor*(rand(3,4)-0.5);
            
        case 2
             r_gt=factor*(rand(3,4) -0.5); 
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
                z= 2*rand -1;
                r=sqrt(1-z^2);
                rigMove=rigLength*[r*cos(phi) ; r*sin(phi); z];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            s_gt=factor*(rand(3,5)-0.5);
            
        case 3
            %3 - 6 mics, 2 sounds, 2D Mics - 3D Sounds, calibrated 
             r_gt=factor*(rand(2,6) -0.5); 
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
             
                rigMove=rigLength*[cos(phi) ; sin(phi); ];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            r_gt=[r_gt; zeros(1,6)];
            
            s_gt=factor*(rand(3,2)-0.5);
            
        case 4
            %4 - 8 mics, 2 sounds, 2D Mics - 3D Sounds, UNcalibrated 
                       
             r_gt=factor*(rand(2,8) -0.5); 
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
             
                rigMove=rigLength*[cos(phi) ; sin(phi); ];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            
              r_gt=[r_gt; zeros(1,8)];
            s_gt=factor*(rand(3,2)-0.5);
            
        case 5
            %5 - 4 mics, 3 sounds, 2D Mics - 3D Sounds, calibrated
            r_gt=factor*(rand(2,4) -0.5);
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
                
                rigMove=rigLength*[cos(phi) ; sin(phi); ];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            s_gt=factor*(rand(3,3)-0.5);
            
        case 6
            %6 - 6 mics, 3 sounds, 3D Mics - 2D Sounds, calibrated
             r_gt=factor*(rand(3,6) -0.5); 
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
                z= 2*rand -1;
                r=sqrt(1-z^2);
                rigMove=rigLength*[r*cos(phi) ; r*sin(phi); z];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            s_gt=factor*(rand(2,3)-0.5);
            
        case 7
            %7 - 8 mics, 3 sounds, 3D Mics - 2D Sounds, UNcalibrated
           
            r_gt=factor*(rand(3,8) -0.5);
            
            %creating receivers on a rig, fixing mic 2, 4 ,6 etc to fit
            for kk=1:size(r_gt,2)/2
                phi=2*pi*rand;
                z= 2*rand -1;
                r=sqrt(1-z^2);
                rigMove=rigLength*[r*cos(phi) ; r*sin(phi); z];
                r_gt(:,kk*2) = r_gt(:,kk*2-1) + rigMove;
            end
            
            s_gt=factor*(rand(2,3)-0.5);
            
        otherwise
            error('Wrong eperiment type. Please check experiments to run in vector "constats"')
    end
    
    D=zeros(size(r_gt,2),size(s_gt,2));
    for ii=1:size(r_gt,2)
        for jj=1:size(s_gt,2)
            D(ii,jj)=norm(r_gt(:,ii)-s_gt(:,jj));
        end
    end
    
    D=D+normrnd(0,stdErr,size(D));
end