function [CI] = Function_CI_calculator_f05(XYY,X1,X2,method)
%% This is used for calculating CI based on max-height method and area method
% 2022/12/19 By Zijiang Yang
% 2022/12/22 Using max-height method for peak calculation: peak height @X = max peak heiht in X ± δx to make method more robust
% 2023/06/20 Remove baseline correction within the function

% 2023/06/22a Combine "global baseline approach" and "local baseline approach"
%             Local baseline approach only work for area method! (=´∀｀)人(´∀｀=)イェーイ
% 2023/06/22b Allow XY with Y as multiple columns...using for loop...o(╥﹏╥)o

% Input
% XYY: X = wavenumber, YY = signal intensities corresponding to X
% X1, X2 = if a wavenumber1 and wavenumber2, then using area method; if wavenumber and 5, then using peak height method
% e.g.., X1 = [1650 1850] or X1 = [1650 5]
% method: baseline method. if = G, then global baseline correction, if = L, then local baseline correction (Walfridson and Kuttainen, 2022)

%% Preparation

CI = NaN(width(XYY) - 1,1);

%% Seperating
for i = 1:1:width(XYY) - 1

    XY = XYY(:,[1,i+1]);
    %% Processing
    % check order
    if XY(1,1) > XY(1:end)
        flip(XY);end   
    % Prepare XY data
    X      = XY(:,1);
    Y      = XY(:,2:end);
    
    if method == 'G'
        %% CI calculation - global approach
        if X1(2) > X1(1) && X2(2) > X2(1) % then it is a range, not radiu, thus area method
            
            IndexX11  = find(X == min(X1)); % the smaller X1 range
            IndexX12  = find(X == max(X1)); % the larger  X1 range
            XY1       = XY(IndexX11:IndexX12,:);
            A1        = trapz(XY1(:,1),XY1(:,2:end),1);
            IndexX21  = find(X == min(X2)); % the smaller X2 range
            IndexX22  = find(X == max(X2)); % the larger  X2 range
            XY2       = XY(IndexX21:IndexX22,:);
            A2        = trapz(XY2(:,1),XY2(:,2:end),1);
        
            CI(i,:)   = (A1 ./ A2)';
        
        elseif X1(2) < X1(1) && X2(2) < X2(1) % radiu is always smaller than wavenumber, thus, max height method
        
            x1 = [X1(1)-X1(2) X1(1)+X1(2)];
            x2 = [X2(1)-X2(2) X2(1)+X2(2)];
        
            IndexX11  = find(X == min(x1)); % the smaller X1 range
            IndexX12  = find(X == max(x1)); % the larger  X1 range
            XY1       = XY(IndexX11:IndexX12,:);
            H1        = max(XY1(:,2:end));
            IndexX21  = find(X == min(x2)); % the smaller X2 range
            IndexX22  = find(X == max(x2)); % the larger  X2 range
            XY2       = XY(IndexX21:IndexX22,:);
            H2        = max(XY2(:,2:end));
        
            CI(i,:)   = (H1 ./ H2)';
        
        else        
            msgbox('Something wrong with X1 and X2, check X1 and X2')            
        end    
    
    elseif method == 'L'
    
        %% CI calculation 
        if X1(2) > X1(1) && X2(2) > X2(1) % then it is a range, not radius, thus area method
    
            % Define base lines
            base1 = linspace(Y(X==X1(1)), Y(X==X1(2)), find(X==X1(2))-find(X==X1(1))+1)';
            base2 = linspace(Y(X==X2(1)), Y(X==X2(2)), find(X==X2(2))-find(X==X2(1))+1)';
        
            % Calculate areas
            A1 = trapz(X(X>=X1(1) & X<=X1(2)), Y(X>=X1(1) & X<=X1(2)) - base1, 1);
            A2 = trapz(X(X>=X2(1) & X<=X2(2)), Y(X>=X2(1) & X<=X2(2)) - base2, 1);
        
            % Calculate contrast index
            CI(i,:) = (A1 ./ A2)';
        
            % Plot for verification
            %     figure
            %     plot(X, Y, 'k');
            %     hold on;
            %     fill([X(X>=X1(1) & X<=X1(2)); flip(X(X>=X1(1) & X<=X1(2)))], [Y(X>=X1(1) & X<=X1(2)); flip(base1)], 'r', 'FaceAlpha', 0.3);
            %     fill([X(X>=X2(1) & X<=X2(2)); flip(X(X>=X2(1) & X<=X2(2)))], [Y(X>=X2(1) & X<=X2(2)); flip(base2)], 'b', 'FaceAlpha', 0.3);
            %     hold off;
    
        elseif X1(2) < X1(1) && X2(2) < X2(1) % radius is always smaller than wavenumber, thus, max height method
    
               CI(i,:)   = NaN;
        else    
        msgbox('Something wrong with X1 and X2, check X1 and X2')        
        end
    end
end
end
