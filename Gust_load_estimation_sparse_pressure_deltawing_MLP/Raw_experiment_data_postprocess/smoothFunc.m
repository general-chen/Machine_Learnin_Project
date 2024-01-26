function [Fx_smoothed] = smoothFunc(Fx,input2)
%    Add some smooth funk to your signals...
% Smooths each column of an input variable via either Savitzky-Golay 
% smoothing or a strict low-pass filter.
%  * This can handle 2D data by evaluating each column at a time.
%  * If input2 is numeric, it defines % the S-G. window size.
%  * If input2 is a string, it is the cutoff frequency of the filter.
%  * If there is no 2nd input, assume 200 pt. S-G smoothing.

[rows,columns] = size(Fx);
Fx_smoothed = zeros(rows,columns);

% Default sgolay; default 200 point window
if nargin == 1
    window = 200;
    for ic = 1:columns
        Fx_smoothed(:,ic)=smooth(Fx(:,ic),window,'sgolay',4);
    end

elseif nargin == 2
    % Default sgolay; input 'window'
    if isnumeric(input2)
        window = input2;
        for ic = 1:columns
            Fx_smoothed(:,ic)=smooth(Fx(:,ic),window,'sgolay',4);
        end
        
    % Other filter:
    else
        FF = str2double(input2);
        % input
        Fs = 1000;              % Sampling Frequency
        Fpass = 0;             % Passband Frequency
%         Fstop = 10;             % Stopband Frequency
        Fstop = FF;             % Stopband Frequency
        Dpass = 0.057501127785; % Passband Ripple
        Dstop = 0.0001;         % Stopband Attenuation
        flag = 'noscale';       % Sampling Flag

        % Calculate the order from the parameters using KAISERORD.
        [N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(Fs/2),[1 0],[Dstop Dpass]);

        %Calculate the coefficients using the FIR1 function.
        num = fir1(N,Wn,TYPE, kaiser(N+1,BETA), flag); % Filter parameters
        GD = N/2;               % Offset delay created by filtering
        clear BETA Dpass Dstop Fpass Fs Fstop TYPE Wn

        % GD used to reset the filtered data to the same time as the input 
        
        for ic = 1:columns
            y = [Fx(:,ic); Fx(end,ic)*ones(ceil(GD),1)];   % input data
            yf = filter(num,1,y);  % filtered output data
            Fx_smoothed(:,ic) = yf(ceil(GD)+1:end); % correct filter shift
        end
        
%         figure; freqz(num)
        

    end
    
else
    error(' Too many input arguments')
end

end

