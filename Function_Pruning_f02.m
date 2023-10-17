function [XYd0,XYd1,XYd2] = Function_Pruning_f02(XY,dv,sp,ck,qt,ws,ss,rg,ct,sm,od)
%% This is used for processing spectrum data under a standardized procesure
% 2022/11/17 Created by Zijiang Yang

%% Parameters
% dv  = 0.1;  % accuracy of spectrum 1D-interpolation
% sp  = 11;   % smoothing range
% ck  = 0;    % check baseline correaction
% qt  = 0;    % Quantile for localized baseline; Lenz et al., 2015 qt = 0
% ws  = 130;  % Window size for smoothing; Lenz et al., 2015 ws = 128
% ss  = 130;  % Step size for smoothing;   Lenz et al., 2015 ss = 128
% 
% rg  = [400 4000]; % range of measured spectrum
% ct  = [650   700; 
%        2250  2450; 
%        000   500; 
%        3500  8000]; % For FTIR [650 700; 2250 2450; 000 500; 3500 8000], Nakano et al., 2021
% 
%% Options
% sm  = 1;    % sm = 0 for moving average (Nakano et al., 2021); sm = 1 for Savitzky-Golay (Lenz et al., 2015)
% od  = 2;    % polynomial order for Savitzky-Golay

%% Calculations
% (1) order check
if XY(1,1) > XY(end,1)
   XY = flip(XY); end

% (2) Alignment of X 
% 1D-interpolation preserved the information from original spectrum
Xvq  = (round(XY(1,1)):dv:round(XY(end,1)))';
Yvq  = interp1(XY(:,1),XY(:,2:end),Xvq,'makima'); % Akima piecewise interpolation, Akima, 1970.
XYvq = [Xvq Yvq];

% (3) Whole range check
% Remove extra range, i.e., we want 400-4000
% Xrg  = (rg(1):dv:rg(2))';
index1 = find(XYvq(:,1) == rg(1)); 
index2 = find(XYvq(:,1) == rg(2));
XYrg   = XYvq(index1:1:index2,:);

% (4) Smoothing
% Set two options: Savitzky-Golay or moving average. define length
if sm == 0
    Ysm = movmean(XYrg(:,2:end),sp); 
elseif sm == 1
    Ysm = sgolayfilt(XYrg(:,2:end),od,sp);
end

XYsm = [XYrg(:,1) Ysm];

% (5) Cutting - remove range that is not used for identification
indexsmp  = NaN([length(ct),2]);
indexcut  = [];


for i = 1:1:height(ct)

    if ct(i,2) < min(XYsm(:,1)) || ct(i,2) > max(XYsm(:,1))
        index(i,1) = NaN;
    else

        if ct(i,1) < min(XYsm(:,1))
            indexsmp(i,1) = find(XYsm(:,1) == min(XYsm(:,1)));
        else
            indexsmp(i,1) = find(XYsm(:,1) == ct(i,1));
        end
    
        if ct(i,2) > max(XYsm(:,1))
            indexsmp(i,2) = find(XYsm(:,1) == max(XYsm(:,1)));
        else
            indexsmp(i,2) = find(XYsm(:,1) == ct(i,2));
        end
    
        indexcut = [indexcut; (indexsmp(i,1):1:indexsmp(i,2))'];

    end
end

XYct  = XYsm;
XYct(indexcut,:) = [];

% (6) Minimal correction
% Let the min spectrum signal as zero
XYmc = XYct;
XYmc(:,2:end) =  XYct(:,2:end) - min(XYct(:,2:end));

% (7) Baseline correction
Ybc   = msbackadj(XYmc(:,1),XYmc(:,2:end),'WindowSize',ws,'StepSize',ss,'EstimationMethod','quantile','RegressionMethod','pchip','QuantileValue',qt,'ShowPlot',ck); 
if ck > 0 
    set(gca, 'XDir','reverse');
end
XYbc = [XYmc(:,1) Ybc];

% (8) Re-Smoothing for cutted range - to prevent extreme values when take diff.
Xbcvq  = ((XYbc(1,1)):dv:(XYbc(end,1)))';
Ybcvq  = interp1(XYbc(:,1),XYbc(:,2:end),Xbcvq,'makima'); % Akima piecewise interpolation, Akima, 1970.
XYsm2  = [Xbcvq Ybcvq];

% (9) Spectrum transformation
XYsm2d0= XYsm2;
Ysm2d0 = XYsm2d0(:,2:end);
% First order derivative
Ysm2d1 = diff(Ysm2d0,1,1);
XYsm2d1= [XYsm2d0(1:end-1,1) Ysm2d1];
XYsm2d1= [XYsm2d1; XYsm2d1(end,:)];
XYsm2d1(:,1)= XYsm2d0(:,1);

% Second order derivative
Ysm2d2 = diff(Ysm2d1,1,1);
XYsm2d2= [XYsm2d1(1:end-2,1) Ysm2d2];
XYsm2d2= [XYsm2d2; XYsm2d2(end,:); XYsm2d2(end,:)];
XYsm2d2(:,1)= XYsm2d0(:,1);

% (10) Re-cutting
X      = XYsm2d0(:,1);
indexsmp  = NaN([length(ct),2]);
indexcut  = [];

for i = 1:1:height(ct)


    if ct(i,1) > max(X) || ct(i,2) < min(X)
        indexsmp(i,1) = NaN;
        indexsmp(i,2) = NaN;
        indexcut = indexcut;
    else
        if ct(i,1) < min(X)
            indexsmp(i,1) = find(X == min(X));
        else
            indexsmp(i,1) = find(X == ct(i,1));
        end

        if ct(i,2) > max(X)
            indexsmp(i,2) = find(X == max(X));
        else
            indexsmp(i,2) = find(X == ct(i,2));
        end

        indexcut = [indexcut; (indexsmp(i,1):1:indexsmp(i,2))'];
    end


    % indexcut = [indexcut; (indexsmp(i,1):1:indexsmp(i,2))'];

end

XYd0  = XYsm2d0;
XYd0(indexcut,:) = [];

XYd1  = XYsm2d1;
XYd1(indexcut,:) = [];

XYd2  = XYsm2d2;
XYd2(indexcut,:) = [];

end




















