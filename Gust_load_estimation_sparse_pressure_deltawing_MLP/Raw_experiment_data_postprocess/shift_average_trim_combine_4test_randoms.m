%% phase average Individual Runs
%  first shift, then phase average, then trim, last combine paressure

clc
clear all
close all

set(0,'DefaultFigureWindowStyle','normal');
startCase = 1;
totalCase = 32;
nRuns     = 10;
data       = cell(totalCase,1);

cross_validation = 8; % number of the folders for cross validation


%% shift, use 1st run as reference to shift
for caseNumber = startCase:totalCase

    % load data
    caseNo     = caseNumber;         % choose case number
    casename   = ['case', num2str(caseNo,'%02i'), '_IdvRuns.mat'];

    data{caseNumber}       = load(casename);

    % shift, use 1st run as reference to shift
    F_Filtered_shift{caseNumber,1}(:,:)      = data{caseNumber}.F_Filtered{1}(:,:);
    F_sStar_shift{caseNumber,1}              = data{caseNumber}.F_sStar{1};

    P_1_8_Filtered_shift{caseNumber,1}(:,:)  = data{caseNumber}.P_1_8_Filtered{1}(:,:);
    P_1_8_sStar_shift{caseNumber,1}          = data{caseNumber}.P_1_8_sStar{1};

    P_9_16_Filtered_shift{caseNumber,1}(:,:) = data{caseNumber}.P_9_16_Filtered{1}(:,:);
    P_9_16_sStar_shift{caseNumber,1}         = data{caseNumber}.P_9_16_sStar{1};

    for noRun =2:nRuns    % different runs; use 1st run as reference to shift
        % shift for force cl cd cm
        for noCol = 1:3   

            [c,lags] = xcorr(data{caseNumber}.F_Filtered{1}(:,noCol),data{caseNumber}.F_Filtered{noRun}(:,noCol));
            [c_max, c_max_index] = max(c);
            shift(noCol) = lags(c_max_index);
            F_Filtered_shift{caseNumber,noRun}(:,noCol) = circshift(data{caseNumber}.F_Filtered{noRun}(:,noCol),shift(noCol));

        end
        % shift for port 1-8
        for noCol = 1:8
            [c,lags] = xcorr(data{caseNumber}.P_1_8_Filtered{1}(:,noCol),data{caseNumber}.P_1_8_Filtered{noRun}(:,noCol));
            [c_max, c_max_index] = max(c);
            shift(noCol) = lags(c_max_index);
            P_1_8_Filtered_shift{caseNumber,noRun}(:,noCol) = circshift(data{caseNumber}.P_1_8_Filtered{noRun}(:,noCol),shift(noCol));
        end
        % shift for port 9-16
        for noCol = 1:8
            [c,lags] = xcorr(data{caseNumber}.P_9_16_Filtered{1}(:,noCol),data{caseNumber}.P_9_16_Filtered{noRun}(:,noCol));
            [c_max, c_max_index] = max(c);
            shift(noCol) = lags(c_max_index);
            P_9_16_Filtered_shift{caseNumber,noRun}(:,noCol) = circshift(data{caseNumber}.P_9_16_Filtered{noRun}(:,noCol),shift(noCol));
        end
        % shift for F_sStar
        [c,lags] = xcorr(data{caseNumber}.F_sStar{1},data{caseNumber}.F_sStar{noRun});
        [c_max, c_max_index] = max(c);
        shift_f_sStar = lags(c_max_index);
        F_sStar_shift{caseNumber,noRun} = circshift(data{caseNumber}.F_sStar{noRun},shift_f_sStar);

        % shift for P_1_8_sStar
        [c,lags] = xcorr(data{caseNumber}.P_1_8_sStar{1},data{caseNumber}.P_1_8_sStar{noRun});
        [c_max, c_max_index] = max(c);
        shift_P_1_8 = lags(c_max_index);
        P_1_8_sStar_shift{caseNumber,noRun} = circshift(data{caseNumber}.P_1_8_sStar{noRun},shift_P_1_8);

        % shift for P_9_16_sStar
        [c,lags] = xcorr(data{caseNumber}.P_9_16_sStar{1},data{caseNumber}.P_9_16_sStar{noRun});
        [c_max, c_max_index] = max(c);
        shift_P_9_16 = lags(c_max_index);
        P_9_16_sStar_shift{caseNumber,noRun} = circshift(data{caseNumber}.P_9_16_sStar{noRun},shift_P_9_16);
    end

end


%% phase average
% phase average
for caseNumber = startCase:totalCase
    
    for noRun = 1:nRuns
        numrowsF(noRun)       = size(F_Filtered_shift{caseNumber,noRun},1);
        numrowsF_sStar(noRun) = size(F_sStar_shift{caseNumber,noRun},1);
    
        numrowsP_1_8(noRun) = size(P_1_8_Filtered_shift{caseNumber,noRun},1);
        numrowsP_1_8_sStar_shift(noRun) = size(P_1_8_sStar_shift{caseNumber,noRun},1);
    
        numrowsP_9_16(noRun) = size(P_9_16_Filtered_shift{caseNumber,noRun},1);
        numrowsP_9_16_sStar_shift(noRun) = size(P_9_16_sStar_shift{caseNumber,noRun},1);
    end
    
    F_Sum_Filter = zeros(min(numrowsF),3,nRuns);
    F_sStarSum = zeros(min(numrowsF_sStar),nRuns);
    
    P_1_8_Sum_Filter = zeros(min(numrowsP_1_8),8,nRuns);
    P_1_8_sStarSum = zeros(min(numrowsP_1_8_sStar_shift),nRuns);
    
    P_9_16_Sum_Filter = zeros(min(numrowsP_9_16),8,nRuns);
    P_9_16_sStarSum = zeros(min(numrowsP_9_16_sStar_shift),nRuns);
    
    for noRun = 1:nRuns
        F_sStarSum(:,noRun)      = F_sStar_shift{caseNumber,noRun}(1:min(numrowsF_sStar));
        P_1_8_sStarSum(:,noRun)  = P_1_8_sStar_shift{caseNumber,noRun}(1:min(numrowsP_1_8_sStar_shift));
        P_9_16_sStarSum(:,noRun) = P_9_16_sStar_shift{caseNumber,noRun}(1:min(numrowsP_9_16_sStar_shift));
        for j =1:3
           F_Sum_Filter(:,j,noRun) = F_Filtered_shift{caseNumber,noRun}(1:min(numrowsF),j);
        end
        for jjj = 1:8
           P_1_8_Sum_Filter(:,jjj,noRun)  = P_1_8_Filtered_shift{caseNumber,noRun}(1:min(numrowsP_1_8),jjj);
           P_9_16_Sum_Filter(:,jjj,noRun) = P_9_16_Filtered_shift{caseNumber,noRun}(1:min(numrowsP_9_16),jjj);
        end
    end
    F_Filtered_PhaseAve{caseNumber}   = mean(F_Sum_Filter,3);
    
    P_1_8_Filtered_PhaseAve{caseNumber}    = mean(P_1_8_Sum_Filter,3);
    P_9_16_Filtered_PhaseAve{caseNumber}   = mean(P_9_16_Sum_Filter,3);
    
    F_sS_PhaseAve{caseNumber}      = mean(F_sStarSum,2);
    P_1_8_sS_PhaseAve{caseNumber}  = mean(P_1_8_sStarSum,2);
    P_9_16_sS_PhaseAve{caseNumber} = mean(P_9_16_sStarSum,2);

end

%% %%% check shift & phase average
% for caseNumber = 1:totalCase
%     f1 = figure;
%     f1.Position = [100 100 1800 800];
% 
%     % check f shift
%     subplot(3,1,1)
%     plot(F_Filtered_PhaseAve{caseNumber},'LineWidth',2, 'Color','r')
%     hold on
%     for i = 1:nRuns
%         plot(F_Filtered_shift{caseNumber,i},'LineStyle','--')
%         hold on
%     end
%     title('f shift')
%     xlim([0 12000])
% 
%     % check P_1_8 shift    
%     subplot(3,1,2)
%     plot(P_1_8_Filtered_PhaseAve{caseNumber},'LineWidth',2, 'Color','r')
%     hold on
%     for i = 1:nRuns
%         plot(P_1_8_Filtered_shift{caseNumber,i},'LineStyle','--')
%         hold on
%     end
%     title('p (1-8) shift')
%     xlim([0 12000])
% 
%     % check P_9_16 shift
%     subplot(3,1,3)
%     plot(P_9_16_Filtered_PhaseAve{caseNumber},'LineWidth',2, 'Color','r')
%     hold on
%     for i = 1:nRuns
%         plot(P_9_16_Filtered_shift{caseNumber,i},'LineStyle','--')
%         hold on
%     end
%     title('p (9 16) shift')
%     xlim([0 12000])
% 
%     sgtitle(['case', num2str(caseNumber,'%02i')]) 
%     % save figure
% %     saveas(gcf,['case', num2str(caseNumber,'%02i') '.png'])
% %     close(f1)
% end

%% trim, using trim mtrix from .csv. Here choose 3000 points to trim, s*=10
lenTrim        = 2000;
lenTrim_before = 500-1; % in order to get 50*50=2500 samples
trimAmendMatrix = [0 0 0 2000 2000 2000 0 2000 1000       0 50 0 50 0 0 0 0 0 -110 0          0 -190 -50 -160 -170 -150 -200 -130 -140 -130 -180 -220];
for caseNumber = startCase:totalCase  

    % load the trim matrix from .txt
    casename   = ['case', num2str(caseNumber,'%02i')];
    whichCase  = fopen('trim_index_matrix_by_force.txt');
    caseUsed   = textscan(whichCase, '%s %s','delimiter','|', 'HeaderLines',1);
    fclose(whichCase);
    clear whichCase
    
    caseIndex     = find(string(strtrim(cell2mat(caseUsed{1}))) == string(casename));
    trimParameter = cell2mat(textscan(cell2mat(caseUsed{2}(caseIndex)),'%f %f %f %f %f %f %f'));

    % trim f
    F_trim{caseNumber} = F_Filtered_PhaseAve{caseNumber}(trimParameter(1)-lenTrim_before:trimParameter(1)+lenTrim,:);
    F_sStar_trim{caseNumber} = F_sS_PhaseAve{caseNumber}(trimParameter(1)-lenTrim_before:trimParameter(1)+lenTrim);
    % normalize S using chord length 0.3m
    F_sStar_trim{caseNumber} = (F_sStar_trim{caseNumber} - F_sStar_trim{caseNumber}(1)) / 0.3;

    %trim P_1_8
    P_1_8_trim{caseNumber} = P_1_8_Filtered_PhaseAve{caseNumber}(trimParameter(3)-lenTrim_before:trimParameter(3)+lenTrim,:);
    P_1_8_sStar_trim{caseNumber} = P_1_8_sS_PhaseAve{caseNumber}(trimParameter(3)-lenTrim_before+trimAmendMatrix(caseNumber):trimParameter(3)+lenTrim+trimAmendMatrix(caseNumber));
    % normalize S using chord length 0.3m
    P_1_8_sStar_trim{caseNumber} = (P_1_8_sStar_trim{caseNumber} - P_1_8_sStar_trim{caseNumber}(1)) / 0.3;

    %trim P_9_16
    P_9_16_trim{caseNumber} = P_9_16_Filtered_PhaseAve{caseNumber}(trimParameter(5)-lenTrim_before:trimParameter(5)+lenTrim,:);
    P_9_16_sStar_trim{caseNumber} = P_9_16_sS_PhaseAve{caseNumber}(trimParameter(5)-lenTrim_before:trimParameter(5)+lenTrim);
    % normalize S using chord length 0.3m
    P_9_16_sStar_trim{caseNumber} = (P_9_16_sStar_trim{caseNumber} - P_9_16_sStar_trim{caseNumber}(1)) / 0.3;
end

%% combine pressure % port 11 is the stagnation point, put it to the end of the matrix
for caseNumber = startCase:totalCase
    P_trim{caseNumber}       = [P_1_8_trim{caseNumber} P_9_16_trim{caseNumber}];
    P_sStar_trim{caseNumber} = P_1_8_sStar_trim{caseNumber};

    P_trim_temp              = P_trim{caseNumber}(:,11); 
    P_trim{caseNumber}(:,11:15) = P_trim{caseNumber}(:,12:16);
    P_trim{caseNumber}(:,16)    = P_trim_temp;
end

% % check trim data
% for caseNumber = startCase:totalCase
% %     plot(P_sStar_trim{caseNumber},P_trim{caseNumber})
%     plot(P_trim{caseNumber})
%     hold on
% end


%% %%%%%%%%%%%%%%%%%%%%%%%%% prepare for ML %%%%%%%%%%%%%%%%%%%%%%%%%
% combine the trained and test cases by column
pwd
currentFolder = pwd;
% delete dataset.h5 file before recreate a new one
testSize = 0.2;
allCase  = totalCase;
allTrain = round(totalCase * (1-testSize));

% read case number for different cross-validation folds from .txt

foldFile    = fopen('delta_wing_32cases_split_random.txt');
caseUsed   = textscan(foldFile, '%s %s','delimiter','|', 'HeaderLines',1);
fclose(foldFile);
clear foldFile

for iii = 1:cross_validation

foldname   = ['Fold_', num2str(iii,'%02i')];
caseIndex(iii)  = find(string(strtrim(cell2mat(caseUsed{1}))) == string(foldname));
caseNum{iii} = cell2mat(textscan(cell2mat(caseUsed{2}(caseIndex(iii))),'%f %f %f %f'));

% create random number for test cases and train cases, 
% used for cross validation
aa = linspace(1,totalCase,totalCase);
case_test{iii} = caseNum{iii}; % sort the bb in right order, ignore nn
case_train{iii} = setdiff(aa,case_test{iii}); % subtract 4 test cases, the rest is 28 train cases


dataTrain_F   = F_trim{case_train{iii}(1)};
dataTrain_F_s = F_sStar_trim{case_train{iii}(1)};
dataTrain_P   = P_trim{case_train{iii}(1)};
dataTrain_P_s = P_sStar_trim{case_train{iii}(1)};

dataTest_F   = F_trim{case_test{iii}(1)};
dataTest_F_s = F_sStar_trim{case_test{iii}(1)};
dataTest_P   = P_trim{case_test{iii}(1)};
dataTest_P_s = P_sStar_trim{case_test{iii}(1)};

%commbine train data set
for i = case_train{iii}(2:length(case_train{iii}))
    dataTrain_F   = [dataTrain_F;F_trim{i}];
    dataTrain_F_s = [dataTrain_F_s;F_sStar_trim{i}];
    dataTrain_P   = [dataTrain_P;P_trim{i}];
    dataTrain_P_s = [dataTrain_P_s;P_sStar_trim{i}];
end
%commbine test data set
for i = case_test{iii}(2:length(case_test{iii}))
    dataTest_F   = [dataTest_F;F_trim{i}];
    dataTest_F_s = [dataTest_F_s;F_sStar_trim{i}];
    dataTest_P   = [dataTest_P;P_trim{i}];
    dataTest_P_s = [dataTest_P_s;P_sStar_trim{i}];
end

% save data to .h5 for python
rootDir    = 'D:\work\myExperiment\model\deltaWing\data_exp\dataOrginazed';
cd([rootDir '\postProcess\dataForML'])

savepath = [rootDir '\' 'postProcess' '\' 'dataForML ' '\' 'dataset_4test_fold' num2str(iii,'%02i') '.h5'];

% save trained cases F
h5create(savepath, '/F/train', size(dataTrain_F));
h5write(savepath, '/F/train', dataTrain_F)
% save trained cases s
h5create(savepath, '/F_s/train', size(dataTrain_F_s));
h5write(savepath, '/F_s/train', dataTrain_F_s)
    
% save test cases F
h5create(savepath, '/F/test', size(dataTest_F));
h5write(savepath, '/F/test', dataTest_F)
% save test cases s
h5create(savepath, '/F_s/test', size(dataTest_F_s));
h5write(savepath, '/F_s/test', dataTest_F_s)
    
% save trained cases P
h5create(savepath, '/P/train', size(dataTrain_P));
h5write(savepath, '/P/train', dataTrain_P)
% save trained cases s
h5create(savepath, '/P_s/train', size(dataTrain_P_s));
h5write(savepath, '/P_s/train', dataTrain_P_s)

% save test cases P
h5create(savepath, '/P/test', size(dataTest_P));
h5write(savepath, '/P/test', dataTest_P)
% save test cases s
h5create(savepath, '/P_s/test', size(dataTest_P_s));
h5write(savepath, '/P_s/test', dataTest_P_s)

cd(currentFolder);

pwd

end



% %% for spectrum analysis
% fs = 1000;
% t  = (0.001:1/fs:12.625)';
% F_NoFiltered = data{1,1}.F_NoFiltered{1}(:,1);
% S  = timetable(seconds(t),F_NoFiltered);


% Train_error =  0.0010416443836090022;
% Test_error  =  0.0010633341572968019;
% b = bar([Train_error Test_error],'FaceColor','flat');
% b.CData(2,:) = [0 0.8 0.8];


