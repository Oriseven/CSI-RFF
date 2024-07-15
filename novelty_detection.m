%%
%% =====================================================================================
%%       Filename:  novelty_detection.m 
%%
%%    Description:  novelty detection algorithems: 'iforest','lof','ocsvm','ocsvm','knn'
%%
%%         Author:  Ruiqi Kong 
%%         Email :  <kr020@ie.cuhk.edu.hk>
%%   Organization:  WiNS group @ The chiniese university of hong kong
%%
%%   Copyright (c)  WiNS group @ The chiniese university of hong kong
%% =====================================================================================
%%

function scores_test=novelty_detection(train_xdata,test_xdata,functionname,dis)
    distance_lof=["euclidean","cityblock","chebychev"];
    distance_ocsvm=["mykernel_euclidean","mykernel_manhattan","mykernel_chebyshev","mykernel_euclideanangle","mykernel_hemitianangle"];
    distance_knn=["Euclidean_distance","Manhattan_distance","Chebyshev_distance","Euclidean_angle","Hermitian_angle"];
    scores_test=[];
    if functionname == "iforest"
        % A score value close to 0 indicates a normal observation, and a value close to 1 indicates an anomaly.
        [Mdl,~,~] = iforest(train_xdata); 
        [~,scores_test] = isanomaly(Mdl,test_xdata);
    elseif functionname == "lof"
        if dis<=length(distance_lof)
            [Mdl,~,~] = lof(train_xdata,"Distance",distance_lof(dis));
            [~,scores_test] = isanomaly(Mdl,test_xdata);
        end
    elseif functionname == "ocsvm"
        % Note that a large positive anomaly score indicates an anomaly in ocsvm, whereas a negative score indicates an anomaly in predict of ClassificationSVM.
        % If the class label variable contains only one class, fitcsvm trains a model for one-class classification and returns a ClassificationSVM object. 
        % To identify anomalies, you must first compute anomaly scores by using the resubPredict or predict object function of ClassificationSVM, 
        % and then identify anomalies by finding observations that have negative scores.
        y = ones(size(train_xdata,1),1);
        train_xdata=train_xdata./100;test_xdata=test_xdata./100;
        SVMModel = fitcsvm(train_xdata,y,'KernelFunction',distance_ocsvm(dis));
        [~,scores_test] = predict(SVMModel,test_xdata);
        scores_test = - scores_test;
    elseif functionname == "knn"
        % D = pdist2(X,Y,Distance,'Smallest',K) computes the distance using the metric specified by Distance 
        % and returns the K smallest pairwise distances to observations in X for each observation in Y in ascending order.
        [dist2, ~]=pdist2(train_xdata,test_xdata, str2func(distance_knn(dis)), 'smallest',ceil(sqrt(size(train_xdata,1))), 'sortindices', 1);
        scores_test = mean(dist2).';
    elseif functionname == "dbscan"
        % The two most important parameter values the model takes are 
        % (i) esp, which specifies the distance between two points i.e., how close the data points should be to one another to be considered part of a cluster; 
        % and (ii) min_samples, which specifies the minimum number of neighbors a point should have in a cluster.
        scores_test = zeros(size(test_xdata,1),1);
        for t=1:size(test_xdata,1)
        idx = dbscan([test_xdata(t,:);train_xdata],0.03,5);%,"Distance",str2func(distance)
        % The function also identifies some outliers (an idx value of ?C1 ) in the data.
        scores_test(t) = (idx(1) == -1);
        end
    end

end

%% distance functions
function D=Euclidean_distance(Y,X)
[nx,~]=size(X);
X = X(:,1:52)+1i*X(:,53:end);Y = Y(:,1:52)+1i*Y(:,53:end);
for i = 1
    dsq = zeros(nx,1,'double');
    for n=1:nx
        dsq(n)=norm(X(n,:) -Y(i,:));
    end
    D = dsq;
end
end


function D=Manhattan_distance(Y,X)
[nx,~]=size(X);
X = X(:,1:52)+1i*X(:,53:end);Y = Y(:,1:52)+1i*Y(:,53:end);
for i = 1
    dsq = zeros(nx,1,'double');
    for n=1:nx
        dsq(n)=sum(abs(X(n,:) -Y(i,:)));
    end
    D = dsq;
end
end

function D=Chebyshev_distance(Y,X)
[nx,~]=size(X);
X = X(:,1:52)+1i*X(:,53:end);Y = Y(:,1:52)+1i*Y(:,53:end);
for i = 1
    dsq = zeros(nx,1,'double');
    for n=1:nx
        dsq(n)=max(abs(X(n,:) -Y(i,:)));
    end
    D = dsq;
end
end


function D=Hermitian_angle(Y,X)
[nx,~]=size(X);
X = X(:,1:52)+1i*X(:,53:end);Y = Y(:,1:52)+1i*Y(:,53:end);
for i = 1
    dsq = zeros(nx,1,'double');
    for n=1:nx
    dsq(n)=abs(sum(X(n,:)'.*Y(i,:).'));   
    dsq(n) = 1 - dsq(n) / (norm(X(n,:))*norm(Y(i,:)));
    end
    D = dsq;
end
end

function D=Euclidean_angle(Y,X)
[nx,~]=size(X);
X = X(:,1:52)+1i*X(:,53:end);Y = Y(:,1:52)+1i*Y(:,53:end);
for i = 1
    dsq = zeros(nx,1,'double');
    for n=1:nx
        num =sum(real(X(n,:)).*real(Y(i,:))) + sum(imag(X(n,:)).*imag(Y(i,:)));
        den = norm(X(n,:))*norm(Y(i,:));
        dsq(n) = abs(1 - num / den);
    end
    D = dsq;
end
end
