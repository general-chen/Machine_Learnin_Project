%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output .csv file
% no filter; just ensemble average, then go to python to filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Notice!!!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CXall: force in x direction (chord direction, parallel to chord)
% CYall: force in y direction (spanwise,)
% CNall: force in z direction (normal to the delta wing)
%
% CLall: lift force (absolute coordinate system, vertical)
% CDall: Drag force (absolute coordinate system, horizontal)
% no moment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

%% read sync exp data. No aoa,just plane force
data_sync = load('../DataSynchronized5.mat');

%% phase average data, and then filter;
% trim data around the gust: 4601-7100, 2500 points in total
trim_start = 4601;
trim_end   = 7100;
total_case = 32; % 32 in total
freq_data  = 1000; 
filter = 400;

for i = 1:total_case
    % phase average force
    CL_temp = mean(data_sync.combineddata{i}.CLall,2);
    data_sync.combineddata{i}.CLall_PhaseAve_filtered = CL_temp(trim_start:trim_end,:);

    CD_temp = mean(data_sync.combineddata{i}.CDall,2);
    data_sync.combineddata{i}.CDall_PhaseAve_filtered = CD_temp(trim_start:trim_end,:);
    
    % phase average pressure
    for j = 1:10
        
        CP_temp = mean(data_sync.combineddata{i}.Cpall,3);
        CP_temp = CP_temp(:, [1:10 12:end 11]); % adjust the order, 11 is the pitot tube
        data_sync.combineddata{i}.Cpall_PhaseAve_filtered = CP_temp(trim_start:trim_end,:,:);

    end

end


%% 
gust_case    = cell(32,1);   % 32 cases in total

for i = 1:total_case  % combine (cp, cl, cd)
    gust_case{i} = [data_sync.combineddata{i}.Cpall_PhaseAve_filtered ...
        data_sync.combineddata{i}.CLall_PhaseAve_filtered ... 
        data_sync.combineddata{i}.CDall_PhaseAve_filtered];

    % transform to cell, to add label, prepare for excel
    gust_case{i}    = num2cell(gust_case{i});

    % create row label
    row_header(1:length(gust_case{1}),1)={['case_', num2str(i,'%02i')]};
    gust_case{i} = [row_header gust_case{i}];
end


%% combine
gust_case_all = gust_case{1};

for j = 2: total_case
    gust_case_all = [gust_case_all; gust_case{j}];
end

%% create column label
col_header = {'Cp1_t'};
for k = [2:10, 12:16, 0]
    col_header_temp = {['Cp' num2str(k) '_t']};
    col_header = [col_header col_header_temp];     %Row cell array (for column labels)
end

col_header = ['case_number' col_header, {'CL'}, {'CD'}];

output_matrix = [col_header; gust_case_all];
writematrix(output_matrix, 'gust_dataframe.csv');  %Write data and both headers









