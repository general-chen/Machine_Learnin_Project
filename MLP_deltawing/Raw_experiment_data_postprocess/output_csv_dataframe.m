%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output .csv file
% insert 1st and 2nd order derivative of all 16 CPs
% cp'(t) = (3CP(t) - 4*CP(t-1) + CP(t-2)) / (2*deltaT)
% cp''(t) = (2CP(t) - 5*CP(t-1) + 4CP(t-2) - CP(t-3)) / (deltaT)^2
% create data frame to excel
% add row and column labels:
% CP(t) CL CD CP'(t) CP''(t)
%  16    1  1   16     16
% 16+1+1+16+16= 50 columns
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

%% read data. No aoa, just plane force
data_sync = load('../DataSynchronized5.mat');

%% phase average data, and then filter;
% trim data around the gust: 4601-7100, 2500 points in total
filter = 30;  % 30 / 1000 = 3%
trim_start = 4601-3; % here -3 is used to calculate the derivative, need 3 ghost points
trim_end   = 7100;
total_case = 32; % 32 in total
freq_data  = 1000; 

for i = 1:total_case
    
    CL_temp = smoothFunc(mean(data_sync.combineddata{i}.CLall,2), num2str(filter));
    data_sync.combineddata{i}.CLall_PhaseAve_filtered = CL_temp(trim_start:trim_end,:);

    CD_temp = smoothFunc(mean(data_sync.combineddata{i}.CDall,2), num2str(filter));
    data_sync.combineddata{i}.CDall_PhaseAve_filtered = CD_temp(trim_start:trim_end,:);
    
    % phase average pressure
    for j = 1:10
        
        CP_temp = smoothFunc(mean(data_sync.combineddata{i}.Cpall,3), num2str(filter));
        CP_temp = CP_temp(:, [1:10 12:end 11]); % adjust the order, 11 is the pitot tube
        data_sync.combineddata{i}.Cpall_PhaseAve_filtered = CP_temp(trim_start:trim_end,:,:);

    end

end


%% calculate 1st and 2nd order derivative
CP_1st_order = cell(32,1); % 1st order derivative of 32 cases
CP_2nd_order = cell(32,1); % 2nd order derivative of 32 cases
gust_case    = cell(32,1);   % 32 cases in total

for i = 1:total_case  % combine (cp, cl, cd)
    gust_case{i} = [data_sync.combineddata{i}.Cpall_PhaseAve_filtered ...
        data_Frieder.combineddata{i}.CLall_PhaseAve_filtered ... 
        data_Frieder.combineddata{i}.CDall_PhaseAve_filtered];

    % insert 1st order derivative here 
    % cp'(t) = (3CP(t) - 4*CP(t-1) + CP(t-2)) / (2*deltaT)

    for ii = 1:16 % calculate derivative of all 16 CPs    
        CP_1st_order{i}(1,ii) = ( 3*gust_case{i}(1,ii) ...
                                - 4*gust_case{i}(1,ii) ...
                                + 1*gust_case{i}(1,ii)) ...
                                / ( 2/freq_data );
        CP_1st_order{i}(2,ii) = ( 3*gust_case{i}(2,ii) ...
                                - 4*gust_case{i}(1,ii) ...
                                + 1*gust_case{i}(1,ii)) ...
                                / ( 2/freq_data );
        for iii = 3:length(gust_case{i})
            CP_1st_order{i}(iii,ii) = ( 3*gust_case{i}(iii,ii) ...
                                      - 4*gust_case{i}(iii-1,ii) ...
                                      + 1*gust_case{i}(iii-2,ii)) ...
                                      / ( 2/freq_data );
        end
    end

    % delete the first 3 ghost points; at first, need to delete only 2
    % ghost points, but in order to keep consistent with 2nd order
    % derivative, here simply delete 3 ghost points
    CP_1st_order{i} = CP_1st_order{i}(4:end,:); 

    % insert 2nd order derivative here 
    % cp''(t) = (2CP(t) - 5*CP(t-1) + 4CP(t-2) - CP(t-3)) / (deltaT)^2

    for jj = 1:16 % calculate derivative of all 16 CPs    
        CP_2nd_order{i}(1,jj) = ( 2*gust_case{i}(1,jj) ...
                                - 5*gust_case{i}(1,jj) ...
                                + 4*gust_case{i}(1,jj) ...
                                - 1*gust_case{i}(1,jj)) ...
                                / ( 1/freq_data )^2;
        CP_2nd_order{i}(2,jj) = ( 2*gust_case{i}(2,jj) ...
                                - 5*gust_case{i}(1,jj) ...
                                + 4*gust_case{i}(1,jj) ...
                                - 1*gust_case{i}(1,jj)) ...
                                / ( 1/freq_data )^2;
        CP_2nd_order{i}(3,jj) = ( 2*gust_case{i}(3,jj) ...
                                - 5*gust_case{i}(2,jj) ...
                                + 4*gust_case{i}(1,jj) ...
                                - 1*gust_case{i}(1,jj)) ...
                                / ( 1/freq_data )^2;
        for jjj = 4:length(gust_case{i})
            CP_2nd_order{i}(jjj,jj) = ( 2*gust_case{i}(jjj,jj) ...
                                      - 5*gust_case{i}(jjj-1,jj) ...
                                      + 4*gust_case{i}(jjj-2,jj) ...
                                      - 1*gust_case{i}(jjj-3,jj)) ...
                                      / ( 1/freq_data )^2;
        end
    end

    % delete the first 3 ghost points
    CP_2nd_order{i} = CP_2nd_order{i}(4:end,:); 

    % delete the first 3 ghost points, to get 2500 points
    gust_case{i} = gust_case{i}(4:end,:);

    % transform to cell, to add label, prepare for excel
    gust_case{i}    = num2cell(gust_case{i});
    CP_2nd_order{i} = num2cell(CP_2nd_order{i});
    CP_1st_order{i} = num2cell(CP_1st_order{i});

    % create row label
    row_header(1:length(gust_case{1}),1)={['case_', num2str(i,'%02i')]};
    gust_case{i} = [row_header gust_case{i}];
end


%% combine
gust_case_all    = gust_case{1};
CP_1st_order_all = CP_1st_order{1};
CP_2nd_order_all = CP_2nd_order{1};
for j = 2: total_case
    gust_case_all    = [gust_case_all; gust_case{j}];
    CP_1st_order_all = [CP_1st_order_all; CP_1st_order{j}];
    CP_2nd_order_all = [CP_2nd_order_all; CP_2nd_order{j}];
end

%% create column label
col_header = {'Cp1_t'};
for k = [2:10, 12:16, 0]
    col_header_temp = {['Cp' num2str(k) '_t']};
    col_header = [col_header col_header_temp];     %Row cell array (for column labels)
end

% create column label for 1st order derivative
col_header_1st = {'Cp1_1st_order_t'};
for k = [2:10, 12:16, 0]
    col_header_temp = {['Cp' num2str(k) '_1st_order_t']};
    col_header_1st = [col_header_1st col_header_temp];     %Row cell array (for column labels)
end

% create column label for 2nd order derivative
col_header_2nd = {'Cp1_2nd_order_t'};
for k = [2:10, 12:16, 0]
    col_header_temp = {['Cp' num2str(k) '_2nd_order_t']};
    col_header_2nd = [col_header_2nd col_header_temp];     %Row cell array (for column labels)
end

col_header = ['case_number' col_header, {'CL'}, {'CD'}, col_header_1st, col_header_2nd];

output_matrix = [col_header; gust_case_all CP_1st_order_all CP_2nd_order_all];

writematrix(output_matrix, 'gust_dataframe.csv');  %Write data and both headers








