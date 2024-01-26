%% Processing Individual Runs

% This script reads the Force and Pressure data collected from Gust Onset
% experiments. It synchronizes the position, force, and pressure data for
% each run of the input test case and saves the compiled results. 
% Phase-averaging is NOT done in this script.
%
%   Prerequisite scripts: smoothFunc;
%   Required files: data folders; and 'amendMatrix.txt'
%
% case name:
% constant : < case01 case02 case03 case04 case05 case06 case07 case08 case09 >
% rampUp   : < case10 case11 case12 case13 case14 case15 case16 case17 case18...
%              case19 case20 >
% rampDown : < case21 case22 case23 case24 case25 case26 case27 case28 case29...
%              case30 case31 case32>
%
% need to amend the angle of attack using the amend matrix.

%% begin

clc
clear all
close all

%% parameter
Nruns = 10; % number of runs
pCalib = 35; % calibration pressure is 35pa
nports       = 8; % number of pressure ports
sensorRangeP = 6895; % the new pressure sensor range is 1 psi = 6895pa.
sensorRangeV = 5; % 0 - 5 voltage

case_all = 32; % number of cases in total
aoa      = cell(case_all,Nruns,1);

for caseNumber = 1:case_all
        close all
        clearvars -except caseNumber aoa Nruns pCalib nports sensorRangeP sensorRangeV
    
    set(0,'DefaultFigureWindowStyle','docked');
    
    caseNo     = caseNumber;         % choose case number
    
    %%%%%%%%%%%%%%%%%%%%%%     O P T I O N S      %%%%%%%%%%%%%%%%%%%%%%%%%
    
    Normalized = true;      % Normalize force and pressure? (default true)
    Filtered   = true;      % Apply filtering to output?    (default false)
    SaveIdv    = true;      % Save individual run data?
    PhaseAve   = true;      % phase average the data?
    Calibrated = true;      % calibrate the pressure data using 35pa?
                            % caused by the traverse
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% case to choose
    
    % Where is the data folder? Leave empty if it's in same directory as this script
    rootDir = 'D:\work\myExperiment\model\deltaWing\data_exp\dataOrginazed'; 
    addpath('D:\work\myExperiment\model\deltaWing\data_exp\dataOrginazed')
    % add this path to matlab path, as smooth function is in it
    
    casename   = ['case', num2str(caseNo,'%02i')];
    casefolder = [rootDir '\' casename];
    
    % load the parameters from txt
    runFile    = fopen([rootDir '\' 'amendMatrix.txt']);
    runsUsed   = textscan(runFile, '%s %s','delimiter','|', 'HeaderLines',12);
    fclose(runFile);
    clear runFile
    
    caseIndex     = find(string(strtrim(cell2mat(runsUsed{1}))) == string(casename));
    runsParameter = cell2mat(textscan(cell2mat(runsUsed{2}(caseIndex)),'%f %f %f %f %f'));
    % runsParameter( aoa_init    aoa_added    s_actuator    w0(mm/s)    w0_amend(mm/s) )
    
    
    % Test case:        Check test matrix
    aoa_init   = runsParameter(1);          
    aoa_end    = runsParameter(1) + runsParameter(2);
    T          = 0.3;         % s perion of 1-cos gust
    
    w0         = runsParameter(4);       % mm/s see test matrix
    w0_amend   = runsParameter(5);      % amend the angle of attack
    
    s_actuator = runsParameter(3);        % mm
    aoa_error  = 0;           % error of aoa 1.75
    
    
    l          = 120;         % distance between two actuators
    gust       = 1;           % [ 1  2  3]  just let gust > 0
    trigger    = 5.14;        % m
    endPos     = 11.2;          % end position of the case, m
    %% DAQ parameters
    U = 1.0;        % m/s, constant throughout motion
    c = 0.3;            % m
    rho = 1000;
    g = 9.8;
    Area = c*c;         % m2 where A = 0.5(c)(2c) = c*c
    pDynamic = 0.5 * rho * U^2;
    Fs = 1000;          % Hz
    
    %% read initial data
    
    fprintf(['\n ___ Processing ' casename ' : %i Runs ___ \n'],Nruns)
    
    lightSwitch              = cell(Nruns,1);
    forceRead                = cell(Nruns,2);      % {initial, data}
    force_1_8_Read           = cell(Nruns,2);      % {initial, data}
    force_9_16_Read          = cell(Nruns,2);      % {initial, data}
    pressure_1_8_Read  = cell(Nruns,2);
    pressure_9_16_Read = cell(Nruns,2); % {initial, data} % actually this is the
                                       % data from light sensor to locate the
                                       % start position of gusts
    runStarts_F        = zeros(Nruns,2);     % [run#, serial , date vector]
    runStarts_P_1_8    = zeros(Nruns,2);
    runStarts_P_9_16   = zeros(Nruns,2);
    
    Posn                = cell(Nruns,1);          % {posInterp, posRaw}
    F_sStar             = cell(Nruns,1);
    P_1_8_sStar         = cell(Nruns,1);
    P_9_16_sStar        = cell(Nruns,1);
    Events              = zeros(Nruns,1);
    triggerIndex_P_1_8     = zeros(Nruns,1);
    triggerIndex_p_9_16    = zeros(Nruns,1);

    t                   = linspace(0,0.3,300);
    
    str1_initF = 'Initial_Force_Data_Run_';
    str1_initP = 'Initial_Pressure_Data_Run_';
    str1_F = 'Force_Data_Run_';
    str1_P = 'Pressure_Data_Run_';
    str3 = '.csv';
    str4 = 'force';
    str5 = 'pressure1-8';
    str6 = 'pressure9-16';
    
    fprintf('Importing:')
    
    cd([casefolder '\' str4]);
    
    % %%%%%%%%%%%%%%%%%%%%%%  read force data  %%%%%%%%%%%%%%%%%%%%%%
    for i = 1:Nruns
    
        currentRun = num2str(i,'%02i');
        %fprintf(currentRun)
        %fprintf(' ')
    
        %%%%%% Force data
        %printf('F ')
        nameF1 = [str1_initF,currentRun,str3];
        nameF2 = [str1_F,currentRun,str3];
        forceRead{i,1} = readmatrix(nameF1); % initial force
        forceRead{i,2} = readmatrix(nameF2); % force
        
        %%%%% Pressure data
        %fprintf('P ')
        nameP1 = [str1_P,currentRun,str3];
        lightSwitch{i} = readmatrix(nameP1); %  it is light switch
    
    %%%%%% Position and start time data
        %fprintf('s* ')
    
        fileIDF = fopen(['Force_Data_Run_' currentRun '.csv']);
        fileIDP = fopen(['Pressure_Data_Run_' currentRun '.csv']);
        % Read time data first, located in first row of Force_Data_*.csv
        timeF = textscan(fileIDF, '%f %f %f %f %f %f', 1,'delimiter', ',');
        fclose(fileIDF);
        fclose(fileIDP);
        
        % Get start date/time in serial format
        timeF_serial = datenum(cell2mat(timeF));
        % track all the run starts
        runStarts_F(i,:) = [i timeF_serial];
        
        % Now, read position data, 2nd column of force data
        posRead     = forceRead{i,2}(:,2);
        triggerRead = lightSwitch{i}(:,3);
        
        % Get position data, and interpolate flat points.
        Posn{i} = [];
    
        changeIdx = [1; find(diff(posRead)~=0)]; % find the step changes
        xq = 1:length(posRead);
        posInterp = interp1(changeIdx(1:end-1),posRead(changeIdx(2:end)),...
                    xq,'spline')';
        posInterp = smoothFunc(posInterp,'8');
        Posn{i} = [posInterp, posRead];
    
        % find trigger from light switch %%% DO NOT USE THIS TRIGGER FROM THE SWITCH %%%
        for ijk = 2:length(triggerRead)  %%% IT HAS A DELAY OF 10ms %%%
          if ( triggerRead(ijk) - triggerRead(ijk-1) ) > 3 % can also 6 7 8
            triggerIndex_F(i) = ijk;
            break
          end
        end
    
    
        % Find index of gust trigger event per run
        if gust == 0
            %for steady cases arbitrarily choose 6.00m as the sync point
            triggers = 6.00;
            F_sStar{i} = (posInterp - posInterp(triggerIndex_F(i)))/c;
        else
            F_sStar{i} = posInterp;
        end
    
        % pitch up, find the aoa at every sampling point
        aoa{caseNumber,i}(1:triggerIndex_F(i))     = aoa_init + aoa_error;
        if aoa_init < aoa_end % pitch up
            aoa{caseNumber,i}(triggerIndex_F(i)+1:triggerIndex_F(i)+300) = ...
            atand( ( (l*tand(aoa_init)) + (w0*w0_amend/2)*(t - (T/(2*pi))*sin(2*pi/T*t)) ) / l );
        else  % pitch down
            aoa{caseNumber,i}(triggerIndex_F(i)+1:triggerIndex_F(i)+300) = ...
            atand( tand(aoa_init) - ((w0*w0_amend/2)*(t - (T/(2*pi))*sin(2*pi/T*t)) ) / l );
        end
    
        aoa{caseNumber,i}(triggerIndex_F(i)+301:length(triggerRead)) = aoa_end;
    
        
        
        % Events(i) = find( F_sStar{i} > 0.0 ,1,'first'); %%% WRONG DO NOT USE %%%
        
        plot(aoa{caseNumber,i})
        axis([5500 6000 0 30])
        hold on
            
    end
    
    %fprintf('\n\n  Events: mean %.1f (std %0.1f)\n',mean(Events),std(Events))
    
    % change folder
    cd(rootDir);
    pwd
    cd([casefolder '\' str5]);
    pwd
% end % for output aoa 

    % %%%%%%%%%%%%%%%%%%%%%%  read pressure 1-8 data  %%%%%%%%%%%%%%%%%%%%%%
    for i = 1:Nruns
    
        currentRun = num2str(i,'%02i');
        %fprintf(currentRun)
    
        %%%%%% Force data it is light sensor
        %fprintf('F ')
        nameF1 = [str1_initF,currentRun,str3];
        nameF2 = [str1_F,currentRun,str3];
        force_1_8_Read{i,2} = readmatrix(nameF2); % force
    
        %%%%% Pressure data
        %fprintf('P ')
        nameP1 = [str1_initP,currentRun,str3];
        nameP2 = [str1_P,currentRun,str3];
        pressure_1_8_Read{i,1} = readmatrix(nameP1); % initial pressure
        pressure_1_8_Read{i,2} = readmatrix(nameP2); % pressure
    
        fileIDF = fopen(['Force_Data_Run_' currentRun '.csv']);
        fileIDP = fopen(['Pressure_Data_Run_' currentRun '.csv']);
        % Read time data first, located in first row of Force_Data_*.csv
        timeP = textscan(fileIDP, '%f %f %f %f %f %f', 1,'delimiter', ',');
        fclose(fileIDP);
    
        % Now, read position data, 2nd column of force data
        posRead     =  pressure_1_8_Read{i,2}(:,2);
        triggerRead =  force_1_8_Read{i,2}(:,3);
        
        % Get position data, and interpolate flat points.
        Posn{i} = [];
    
        changeIdx = [1; find(diff(posRead)~=0)]; % find the step changes
        xq = 1:length(posRead);
        posInterp = interp1(changeIdx(1:end-1),posRead(changeIdx(2:end)),...
                    xq,'spline')';
        posInterp = smoothFunc(posInterp,'8');
        Posn{i} = [posInterp, posRead];
    
        for ijk = 2:length(triggerRead)
          if ( triggerRead(ijk) - triggerRead(ijk-1) ) > 3
            triggerIndex_P_1_8(i) = ijk;
            break
          end
        end
    
    
        % Find index of gust trigger event per run
        if gust == 0
            %for steady cases arbitrarily choose 6.00m as the sync point
            triggers = 6.00;
            P_1_8_sStar{i} = (posInterp - posInterp(triggerIndex_P_1_8(i)))/c;
        else
            P_1_8_sStar{i} = posInterp;
        end
            
    end
    
    %fprintf('\n')
    % change folder
    cd(rootDir);
    pwd
    cd([casefolder '\' str6]);
    pwd
    
    % %%%%%%%%%%%%%%%%%%%%%%  read pressure 9-16 data  %%%%%%%%%%%%%%%%%%%%%%
    for i = 1:Nruns
    
        currentRun = num2str(i,'%02i');
        %fprintf(currentRun)
    
        %%%%%% Force data it is light sensor
        %fprintf('F ')
        nameF1 = [str1_initF,currentRun,str3];
        nameF2 = [str1_F,currentRun,str3];
        force_9_16_Read{i,2} = readmatrix(nameF2); % force
    
        %%%%% Pressure data
        %fprintf('P ')
        nameP1 = [str1_initP,currentRun,str3];
        nameP2 = [str1_P,currentRun,str3];
        pressure_9_16_Read{i,1} = readmatrix(nameP1); % initial pressure
        pressure_9_16_Read{i,2} = readmatrix(nameP2); % pressure
    
        fileIDF = fopen(['Force_Data_Run_' currentRun '.csv']);
        fileIDP = fopen(['Pressure_Data_Run_' currentRun '.csv']);
        % Read time data first, located in first row of Force_Data_*.csv
        timeP = textscan(fileIDP, '%f %f %f %f %f %f', 1,'delimiter', ',');
        fclose(fileIDP);
        
        % Now, read position data, 2nd column of force data
        posRead     =  pressure_9_16_Read{i,2}(:,2);
        triggerRead =  force_9_16_Read{i,2}(:,3);
        
        % Get position data, and interpolate flat points.
        Posn{i} = [];
    
        changeIdx = [1; find(diff(posRead)~=0)]; % find the step changes
        xq = 1:length(posRead);
        posInterp = interp1(changeIdx(1:end-1),posRead(changeIdx(2:end)),...
                    xq,'spline')';
        posInterp = smoothFunc(posInterp,'8');
        Posn{i} = [posInterp, posRead];
    
        for ijk = 2:length(triggerRead)
          if ( triggerRead(ijk) - triggerRead(ijk-1) ) > 3
            triggerIndex_P_9_16(i) = ijk;
            break
          end
        end
    
    
        % Find index of gust trigger event per run
        if gust == 0
            %for steady cases arbitrarily choose 6.00m as the sync point
            triggers = 6.00;
            P_9_16_sStar{i} = (posInterp - posInterp(triggerIndex_P_9_16(i)))/c;
        else
            P_9_16_sStar{i} = posInterp;
        end
            
    end
    
    % change folder
    cd(rootDir);
    pwd
    
    %% Process Forces
    
    % %%%%%%%%%%%%%%%%%%%%%%  process force data  %%%%%%%%%%%%%%%%%%%%%%
    %fprintf('\n\n Processing Forces ... \n')
    
    Forces = cell(Nruns,1);
    F_NoFiltered = cell(Nruns,1);
    InitialF = zeros(Nruns,6);
    lenF = zeros(Nruns,1);
    LF = 0;
    
    % RotAoa =    [cosd(aoa)   0   -sind(aoa);
    %              0          1           0;
    %              sind(aoa)   0    cosd(aoa)];
    
    % >> We end up with:
    %    x in direction of motion
    %     +Fx = -Drag Force;   +Mx = +Roll Moment (right tip down)
    %    z aligned with lift
    %     +Fz = +Lift Force;   +Mz = +Yaw Moment
    %    y to satisfy right handed convention
    %     +Fy = +Side Force;   +My = +Pitching Moment (nose up)
    ForceText = {'Drag Force' 'Side Force' 'Lift Force' ...
                 'Roll  Moment' 'Pitching Moment' 'Yaw Moment'};
         
    for jjj = 1:Nruns
      lenF(jjj) = size(forceRead{jjj,2}, 1);
    end
             
    Force_run   = cell(Nruns,1);
    
    for i2 = 1:Nruns
        
        % transform initial data to tank coordinate system
        RotAoa_init =    [cosd(aoa_init)   0   -sind(aoa_init);
                         0                 1                    0;
                         sind(aoa_init)    0   cosd(aoa_init)];
        Force_initial  = mean(forceRead{i2,1}(:,4:end));
        F_init         = (RotAoa_init*Force_initial(:,1:3)')';
        M_init         = Force_initial(:,4:6);
    
        Force_run{i2}  = forceRead{i2,2}(:,4:end);
    
            %Transform run data to tank coordinates
          for ii2 = 1:min(lenF)
            
              RotAoa{ii2} =    [cosd(aoa{caseNumber,i2}(ii2))   0   -sind(aoa{caseNumber,i2}(ii2));
                                0                 1                    0;
                             sind(aoa{caseNumber,i2}(ii2))    0   cosd(aoa{caseNumber,i2}(ii2))];
            
            F_run(ii2,:) = (RotAoa{ii2}*Force_run{i2}(ii2,1:3)')';
            M_run(ii2,:) = Force_run{i2}(ii2,4:6);
          end
    
        Fxyz = F_run - F_init;
        Mxyz = M_run - M_init;
        
        %Move moments to 1/4 chord from 1/2 chord. It's easiest to move the
        %moment in the wing reference frame rather than tank frame.
        Mxyz(:,1) = Mxyz(:,1) ;
        Mxyz(:,2) = -Mxyz(:,2) - Fxyz(:,3)*(c/4);
    %     Mxyz(:,2) = Mxyz(:,2);
        Mxyz(:,3) = Mxyz(:,3) + Fxyz(:,2)*(c/4);
        
        if Normalized
            Forces{i2} = [ Fxyz Mxyz/c ] ./ (0.5*1000*Area*U^2);
        else
            Forces{i2} = [Fxyz            Mxyz];
            %            [-Drag Side Lift Roll Pitch Yaw]        
        end
    
        LF = max(LF,lenF(i2));
    
        % adjust the order
        F_NoFiltered{i2}(:,1) = Forces{i2}(:,3);     % Cl
        F_NoFiltered{i2}(:,2) = -Forces{i2}(:,1);     % Cd
        F_NoFiltered{i2}(:,3) = Forces{i2}(:,5);     % Cm
    
        % filtered sampling frequency
        F_Filtered{i2}(:,1)   = smoothFunc(Forces{i2}(:,3),'32');     % Cl
        F_Filtered{i2}(:,2)   = -smoothFunc(Forces{i2}(:,1),'32');     % Cd
        F_Filtered{i2}(:,3)   = smoothFunc(Forces{i2}(:,5),'32');     % Cm
    
    
    end
    
    hold off
    
    
    %% Process pressure data
    % %%%%%%%%%%%%%%%%%%%%%%  process pressure 1-8 data  %%%%%%%%%%%%%%%%%%%%%%
    
    %fprintf('Processing Pressure 1-8 data ... \n')
    
    realInitP_1_8 = cell(1,Nruns);
    realP_1_8     = cell(1,Nruns);
    lenP_1_8      = zeros(Nruns,1);
    LF            = 0;
    
    for jjj = 1:Nruns
      lenP_1_8(jjj) = size(pressure_1_8_Read{jjj,2}, 1);
    end
    
    for i = 1:Nruns
        
        realInitP_1_8{i} = mean(pressure_1_8_Read{i,1}(:,4:end)) ...
            * sensorRangeP / sensorRangeV;
        
        realP_1_8{i} = pressure_1_8_Read{i,2}(:,4:end) ...
            * sensorRangeP / sensorRangeV;
        
        if Normalized
            pressure_1_8{i} = (realP_1_8{i} - realInitP_1_8{i}) ./ pDynamic;
        else
            pressure_1_8 = realP_1_8;
        end
    
        if Calibrated
            for k = 1: nports
                pressure_1_8{i}(:,k) = (pressure_1_8{i}(:,k) + pCalib./pDynamic) * (-1);
            end
        end
    
    end
    
    % %%%%%%%%%%%%%%%%%%%%%%  process pressure 9-16 data  %%%%%%%%%%%%%%%%%%%%%%
    
    %fprintf('Processing Pressure 9-16 data ... \n')
    
    realInitP_9_16 = cell(1,Nruns);
    realP_9_16     = cell(1,Nruns);
    lenP_9_16      = zeros(Nruns,1);
    LF             = 0;
    
    for jjj = 1:Nruns
      lenP_9_16(jjj) = size(pressure_9_16_Read{jjj,2}, 1);
    end
    
    for i = 1:Nruns
        
        realInitP_9_16{i} = mean(pressure_9_16_Read{i,1}(:,4:end)) ...
            * sensorRangeP / sensorRangeV;
        
        realP_9_16{i} = pressure_9_16_Read{i,2}(:,4:end) ...
            * sensorRangeP / sensorRangeV;
        
        if Normalized      % transform to cp
            pressure_9_16{i} = (realP_9_16{i} - realInitP_9_16{i}) ./ pDynamic;
        else
            pressure_9_16 = realP_9_16;
        end
    
        if Calibrated
            for k = 1: nports
                pressure_9_16{i}(:,k) = (pressure_9_16{i}(:,k) + pCalib./pDynamic) * (-1);
            end
        end
        
    end
    
    % filter pressure data
    for i100 = 1:Nruns
        P_1_8_NoFiltered{i100}  = pressure_1_8{i100};
        P_9_16_NoFiltered{i100} = pressure_9_16{i100};

        P_1_8_Filtered{i100}    = smoothFunc(P_1_8_NoFiltered{i100},'32');
        P_9_16_Filtered{i100}   = smoothFunc(P_9_16_NoFiltered{i100},'32');

    end
  
%     % do not combine here. combine them after phhase avereage
%     %% combine pressure 1-8 and pressure 9-16,
%     %  and adjust the orders: (port 11 is stagnation pressure)
%     %  put port 11 at the end
%     %  1 2 3 4 5 6 7 8 9 10 12 13 14 15 11
%     
%     P_NoFiltered   = cell(Nruns,1);
%     P_sStar        = cell(Nruns,1);
%     triggerIndex_P = zeros(Nruns,1);
%     for i100 = 1:Nruns
%     
%         [m1, n1] = size(pressure_1_8{i100});
%         [m2, n2] = size(pressure_9_16{i100});
%     
%         m = min(m1,m2);
%         % combine
%         P_NoFiltered{i100} = [pressure_1_8{i100}(1:m,:), pressure_9_16{i100}(1:m,:)];
%         P_sStar{i100}      = (P_1_8_sStar{i100}(1:m,:) + P_9_16_sStar{i100}(1:m,:))/2;
%         % adjust order
%         P_NoFiltered_temp           = P_NoFiltered{i100}(:,11); 
%         P_NoFiltered{i100}(:,11:15) = P_NoFiltered{i100}(:,12:16);
%         P_NoFiltered{i100}(:,16)    = P_NoFiltered_temp;
%     
%         % trigger index p
%         triggerIndex_P(i100) = ( triggerIndex_P_1_8(i100) +triggerIndex_P_9_16(i100) ) / 2; 
%     
%         P_Filtered{i100}   = smoothFunc(P_NoFiltered{i100},'8');  
%     
%     end
    
    
    
    
%     %% no need to sync data here    
%     %% %% %%%%%%%%%%%%%%%%%%%% syncing pressure to force %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % fprintf('synchronizing: \n')
%     % 
%     % p_sync_1_8      = cell(Nruns,1);
%     % p_sync_9_16     = cell(Nruns,1);
%     % inip_1_8_sync   = zeros(Nruns,8);
%     % inip_9_16_sync  = zeros(Nruns,8);
%     % delay_1_8       = zeros(Nruns,1);
%     % delay_9_16      = zeros(Nruns,1);
%     % sync_1_8        = zeros(Nruns,1);
%     % sync_9_16       = zeros(Nruns,1);
%     % trim            = [0 0];
%     % 
%     % fprintf(' aligning data ... ')
%     % 
%     % for i = 1:Nruns
%     % 
%     %     % delay in days, using serial time: < F_time - P_time >
%     %     days_1_8      = runStarts_F(i,2) - runStarts_P_1_8(i,2);
%     %     days_9_16     = runStarts_F(i,2) - runStarts_P_9_16(i,2);
%     %     delay_1_8(i)  = ( 24*60*60 )*days_1_8;   %delay in seconds
%     %     delay_9_16(i) = ( 24*60*60 )*days_9_16;  %delay in seconds
%     %     sync_1_8(i)   = round(delay_1_8(i)*Fs);  %delay in timesteps
%     %     sync_9_16(i)  = round(delay_9_16(i)*Fs); %delay in timesteps
%     % 
%     %     % 'sync' will bring the pressure data in line with the force data for
%     %     % the corresponding run.
%     %     shift_1_8 = sync_1_8(i);
%     %     if shift_1_8 > 0
%     %        % means pressrue starts before force data.
%     %        % only keep data after force starts, pad the end
%     %        padded        = repmat(pressure_1_8{i}(end,:), abs(shift_1_8),1);
%     %        p_sync_1_8{i} = [pressure_1_8{i}(1+shift_1_8:end,:); padded];
%     % 
%     %     elseif shift_1_8 < 0
%     %        % means pressurestarts after force data
%     %        % pad front with first value repeated, trim the end
%     %        padded        = repmat(pressure_1_8{i}(1,:),abs(shift_1_8),1);
%     %        p_sync_1_8{i} = [ padded; pressure_1_8{i}(1:end+shift_1_8,:)];
%     %     end
%     % 
%     % end
%     
%     
%     
%     %% do not phase average here, it is not correct.
%     %% %%%%%%%%%%%%%%%%% phase average data  %%%%%%%%%%%%%%%%%%%%%%%%
%     numrowsF = zeros(Nruns,1);
%     numrowsP = zeros(Nruns,1);
%     for i200 = 1:Nruns
%         numrowsF(i200) = size(F_sStar{i200}, 1);
%         numrowsP(i200) = size(P_sStar{i200}, 1);
%     end
%     
%     F_sStarSum = zeros(min(numrowsF),Nruns);
%     P_sStarSum = zeros(min(numrowsP),Nruns);
%     F_Sum_NoFilter = zeros(min(numrowsF),3,Nruns);
%     P_Sum_NoFilter = zeros(min(numrowsP),nports*2,Nruns);
%     F_Sum_Filter = zeros(min(numrowsF),3,Nruns);
%     P_Sum_Filter = zeros(min(numrowsP),nports*2,Nruns);
%     
%     if PhaseAve
%         for i = 1:Nruns
%             F_sStarSum(:,i) = F_sStar{i}(1:min(numrowsF));
%             P_sStarSum(:,i) = P_sStar{i}(1:min(numrowsP));
%             for j =1:3
%                F_Sum_NoFilter(:,j,i) = F_NoFiltered{i}(1:min(numrowsF),j);
%                F_Sum_Filter(:,j,i) = F_Filtered{i}(1:min(numrowsF),j);
%             end
%             for jjj = 1:nports*2
%                P_Sum_NoFilter(:,jjj,i) = P_NoFiltered{i}(1:min(numrowsP),jjj);
%                P_Sum_Filter(:,jjj,i) = P_Filtered{i}(1:min(numrowsP),jjj);
%             end
%         end
%         F_NoFiltered_PhaseAve = mean(F_Sum_NoFilter,3);
%         F_Filtered_PhaseAve   = mean(F_Sum_Filter,3);
%     
%         P_NoFiltered_PhaseAve = mean(P_Sum_NoFilter,3);
%         P_Filtered_PhaseAve   = mean(P_Sum_Filter,3);
%     
%         F_sS_PhaseAve = mean(F_sStarSum,2);
%         P_sS_PhaseAve = mean(P_sStarSum,2);
%     end
    
    
    pwd
    
    %% save data
    
%     savefile = [rootDir '\' 'postProcess' '\' 'individualRuns' '\' casename '_IdvRuns'];
%     save(savefile,'triggerIndex_F' , 'F_sStar', 'F_NoFiltered', 'F_Filtered',...
%         'triggerIndex_P_1_8', 'P_1_8_sStar', 'P_1_8_NoFiltered', 'P_1_8_Filtered',...
%         'triggerIndex_P_9_16', 'P_9_16_sStar', 'P_9_16_NoFiltered', 'P_9_16_Filtered', '-v7.3')
%     
%     
%     
%     disp([casename,' -- Saved.'])

end











