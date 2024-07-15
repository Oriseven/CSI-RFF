%%
%% =====================================================================================
%%       Filename:  DBSCAN_main.m 
%%
%%    Description:  CSI-RFF implementation for DBSCAN algorithm
%%
%%         Author:  Ruiqi Kong 
%%         Email :  <kr020@ie.cuhk.edu.hk>
%%   Organization:  WiNS group @ The chiniese university of hong kong
%%
%%   Copyright (c)  WiNS group @ The chiniese university of hong kong
%% =====================================================================================
%%
%% dataloader
clear;
load("CSI_data.mat");
% -------------------
% Each row represent CSI from one NIC; Each column represent CSI collected
% in one condition.
% NICs_order =["ESP32C1","ESP32C2","ESP32C3","ESP32C4","ESP32C5",...
%     "AX200C1","AX200C2","AC8260C1","AC7260C1",...
%     "AC7265C1","RTL8812BU","AR9271C1","AR9271C2","AR9271C3","AR9271C4"];
% Conditions_order =
% ["RoomA_static","RoomA_static","RoomA_mobile","RoomA_mobile","RoomB_static","RoomB_static","RoomB_mobile","RoomB_mobile"];

%% fingerprint construction
N_csi = 20; % the number of CSI measurements used for fingerprint construction
N_rx = 1:4; % used rx chains
enable_oe = 1; % enable outlier elimination. (Algorithm 1 in CSI-RFF paper)
n_taps = 8; % the number leakaged taps caused by pulse shaping
fingerprints=Fingerprint(N_csi,N_rx,enable_oe,n_taps);
for nic=1:size(CSI,1) 
    get_micro_csi_group(fingerprints,CSI(nic,:));
end
clearvars -except fingerprints;
%% fingerprint normalization
data=struct2cell(fingerprints.devices);
for i=1:length(data)
    for j= 1:length(data{i,1}{1,1})
        data{i,1}{1,1}{1,j}=zscore((data{i,1}{1,1}{1,j}),[],4);
    end
end
clearvars -except fingerprints data;
%% pairwise distances
distance={'Euclidean_distance','Manhattan_distance','Chebyshev_distance','Euclidean_angle','Hermitian_angle'};
test_enviroment = [5,6;7,8].';distance_value={};self_distance={};
for env = 1:size(test_enviroment,2)
    for dis = 2
        for legal = 1: length(data)
            train_xdata=[];train_ylabel=[];
            f=squeeze(cell2mat(data{legal,1}{1,1}(1,[1:4]).'));
            train_xdata=cat(1,train_xdata,f);
            train_xdata=cat(2,real(train_xdata),imag(train_xdata));
            dist=pdist2(train_xdata,train_xdata, str2func(distance{dis}));
            self_distance{dis,env,legal}=dist;
            for test_device = 1: length(data)
                test_xdata=[];test_ylabel=[];
                f=squeeze(cell2mat(data{test_device,1}{1,1}(1,test_enviroment(:,env)).'));
                test_xdata=cat(1,test_xdata,f);
                test_xdata=cat(2,real(test_xdata),imag(test_xdata));
                dist=pdist2(test_xdata,train_xdata, str2func(distance{dis}));
                distance_value{dis,env,test_device,legal}=dist;
            end
        end
    end
end
clearvars -except fingerprints data self_distance distance_value; 
%% dbscan score
distance={'Euclidean_distance','Manhattan_distance','Chebyshev_distance','Euclidean_angle','Hermitian_angle'};
test_enviroment = [5,6;7,8].';tic;steps=50;EPS=zeros(length(distance),size(test_enviroment,2), length(data),steps);
for env = 1:size(test_enviroment,2)   
    for dis = 2
        for legal = 1: length(data)
            distance_train=self_distance{dis,env,legal};   
            [~, maxv] = bounds(distance_value{dis,env,legal,legal}(:));
            EPS(dis,env,legal,:) = linspace(0.001,maxv,steps);
            for test_device = 1: length(data)
            distance_test=distance_value{dis,env,test_device,legal};
            scores_test = zeros(size(distance_test,1),steps);
                for t=1:size(distance_test,1)
                    D=[[0;distance_test(t,:).'],[distance_test(t,:);distance_train]];
                    for eps = 1:steps
                    [idx, corepts] = dbscan(D,EPS(dis,env,legal,eps),min(53,size(distance_train,1)),'Distance','precomputed');
                    % The function also identifies some outliers (an idx value of -1 ) in the data.
                    scores_test(t,eps) = (idx(1) == -1);
                    end
                end
            scores_dbscan{dis,env,test_device,legal}=scores_test;
            end
        end
    end
end
toc;
%% authentication dbscan
distance={'Euclidean_distance','Manhattan_distance','Chebyshev_distance','Euclidean_angle','Hermitian_angle'};
test_enviroment = [5,6;7,8].';
FRR=[];%false reject rate = false alarm rate
ADR=[];%attack detection rate
threshold = 1/100;
ADR_summary=[];steps=50;
for env = 1:size(test_enviroment,2)
    disp(['-------------env:   ', num2str(env), '--------------']);
    for dis = 2
        disp(['-------------distance:   ', distance{dis}, '--------------']);
        ADR_all=[];
        for legal = 1: length(data)
            for i_nic=1: length(data)
                for eps = 1:steps 
                    if i_nic==legal
                        score=scores_dbscan{dis,env,legal,legal}(:,eps);
                        FRR(dis,env,legal,eps)=sum(score)/length(score);
                    else
                        score=scores_dbscan{dis,env,i_nic,legal}(:,eps);
                        far_num = length(score)- sum(score);
                        test_num = length(score);
                        ADR(dis,env,i_nic,legal,eps) = 1-far_num/test_num;
                    end
                end
            end 
            indx = find(FRR(dis,env,legal,:) <= threshold);
            ADR_all(legal,:)=ADR(dis,env,setdiff(1:15,legal),legal,indx(1))*100;
        end
        ADR_summary(dis,env,:)=[mean(ADR_all(:)),max(ADR_all(:)),min(ADR_all(:))];
        disp(['When FRR=0, average ADR= ' num2str(ADR_summary(dis,env,1)) '  Max ADR= ' num2str(ADR_summary(dis,env,2)) ' MinADR= ' num2str(ADR_summary(dis,env,3))])

    end
end

%% distance
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