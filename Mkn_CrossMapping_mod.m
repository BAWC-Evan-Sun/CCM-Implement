function cr = Mkn_CrossMapping_mod(Data, tau, E, rgns)


    N0 = 10;
    
    n_ts = size(Data, 2);
    n_pts = size(Data, 1);
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
    rgns1 = rgns(:,1);
    rgns2 = rgns(:,2);
    for i = rgns1
        for j = rgns2
           
            X = Data(:,i);
            Y = Data(:,j);
           
            cond_idx = orders;
            cond_idx([i,j]) = [];

            Xests = (squeeze(ests_all(cond_idx, i, :)));
            Yests = (squeeze(ests_all(cond_idx, j, :)));
            XestY = (squeeze(ests_all(j, i, :)));
            YestX = (squeeze(ests_all(i, j, :)));
            XT = X(T+N0-1:end);
            YT = Y(T+N0-1:end);

            % rT = ones(length(XT), 1);
            for z = 1:size(Xests, 2)
                cXest = Xests(:,z);
                cYest = Yests(:,z);
               
               
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
   
                Pyy_cond = (abs(corr(YE_XZ, YT)) - abs(corr(YE_Z, YT)))/abs(corr(YE_XZ, YT)); % In corr_out
   
               
                errY_XZ = YE_XZ - YT;
   
               
                errY_Z = YE_Z - YT;
   
               
                uX2Y = (var(errY_Z) - var(errY_XZ))/var(errY_Z); % In ev_out
   
               
                % ev_out(j,i) = uY2X;
                cr(1, z) = uX2Y;
                cr(2, z) = uY2X;
            end
           
        end
    end
end