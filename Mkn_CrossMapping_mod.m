% Implementation of multivariate-KNN-Search-based conditional cCCM
function cr = Mkn_CrossMapping_mod(Data, tau, E, rgns)
% Inputs: Data - Matrix of input time series (More than two columns). The 
%                columns used to calculate cCCM are indicated by the
%                variable "rgns".
%         tau  - Time lag in cCCM algorithm
%         E    - The dimension of the shadow manifold in cCCM algorithm
%         rgns - An 1 by 2 vector. e.g., if rgns = [1,2], this function
%                take the 1st column as X and 2nd column as Y, and
%                calculate the conditional cCCM between X and Y
%                conditioning on all other columns. 
%
% Outputs: cr  - causality ratios, which is the conditional cCCM value.
    
    N0 = 10; % In this function, we set N0 = 10. This could be adjusted at the user's convenience
    
    n_ts = size(Data, 2);
    n_pts = size(Data, 1);

    if n_ts >= n_pts
        error("Wrong Dimension for Data")
    end
    nrgns = size(rgns, 1);

    cr = zeros(2, nrgns);
    crs = zeros(2, nrgns);
    orders = 1:n_ts;

    T = 1 + (E-1)*tau;
    Lt = n_pts - T + 1;

    ests_all = zeros(n_ts, n_ts, Lt);
    for i = 1:n_ts-1
        for j = i+1:n_ts
           
            X = Data(:,i);
            Y = Data(:,j);
           
            [Pxx, Pyy, XestY, YestX] = cCCM(X, Y, tau, E, N0);
            ests_all(i,j,:) = YestX;
            ests_all(j,i,:) = XestY;
            ccms(i,j) = Pyy;
            ccms(j,i) = Pxx;
        end
    end

    for k = 1:size(rgns,1)
        i = rgns(k,1);
        j = rgns(k,2);
        X = Data(:,i);
        Y = Data(:,j);
       
        cond_idx = orders;
        cond_idx([i,j]) = [];

        Xests = (squeeze(ests_all(cond_idx, i, :)));
        Yests = (squeeze(ests_all(cond_idx, j, :)));
        if size(Xests,2) > size(Xests,1)
            Xests = Xests';
        end
        if size(Yests,2) > size(Yests,1)
            Yests = Yests';
        end
        XestY = (squeeze(ests_all(j, i, :)));
        YestX = (squeeze(ests_all(i, j, :)));
        XT = X(T+N0-1:end);
        YT = Y(T+N0-1:end);

        % rT = ones(length(XT), 1);
        %for z = 1:size(Xests, 2)
        cXest = Xests;
        cYest = Yests;
       
       
        Xe_all = [cXest, XestY];
        Xe_all = Xe_all(N0:end,:);
        cXest = cXest(N0:end,:);
        coef_X_YZ = Xe_all\XT;
        XE_YZ = Xe_all*coef_X_YZ;
        coef_X_Z = cXest\XT;
        XE_Z = cXest*coef_X_Z;

        Pxx_cond = (abs(corr(XE_YZ, XT)) - abs(corr(XE_Z, XT)))/abs(corr(XE_YZ, XT)); % In corr_out
        errX_YZ = XE_YZ - XT;
        errX_Z = XE_Z - XT;
        uY2X = (var(errX_Z) - var(errX_YZ))/var(errX_Z); % In ev_out
        % ev_out(i,j) = uX2Y;

        Ye_all = [cYest, YestX];
        Ye_all = Ye_all(N0:end,:);

        cYest = cYest(N0:end,:);
        coef_Y_XZ = Ye_all\YT;
        YE_XZ = Ye_all*coef_Y_XZ;
        coef_Y_Z = cYest\YT;
        YE_Z = cYest*coef_Y_Z;

        
       
        errY_XZ = YE_XZ - YT;

       
        errY_Z = YE_Z - YT;

       
        uX2Y = (var(errY_Z) - var(errY_XZ))/var(errY_Z); % In ev_out

       
        % ev_out(j,i) = uY2X;
        cr(1, k) = uX2Y;
        cr(2, k) = uY2X;
        %end
           
        
    end
end