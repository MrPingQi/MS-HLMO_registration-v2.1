%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Gao Chenzhong
% Contact: gao-pingqi@qq.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [cor1,cor2] = Multiscale_Matching(descriptors_1,descriptors_2,...
    nOctaves1,nOctaves2,nLayers,Error,scl_flag,par_flag,K)
%% Matching
matches = cell(nOctaves1,nLayers,nOctaves2,nLayers);
confidence = zeros(nOctaves1,nLayers,nOctaves2,nLayers);
if par_flag
    if scl_flag
        nScales = nOctaves1*nLayers*nOctaves2*nLayers;
        k = 1:nScales;
        Octave2 = ceil(k/(nOctaves2*nLayers*nLayers));
        Octave1 = mod(ceil(k/(nLayers*nLayers))-1,nOctaves1)+1;
    else
        nScales = min(nOctaves1,nOctaves2)*nLayers*nLayers;
        k = 1:nScales;
        Octave2 = ceil(k/(nLayers*nLayers));
        Octave1 = Octave2;
    end
        Layer2 = mod(ceil(k/nLayers)-1,nLayers)+1;
        Layer1 = mod(k-1,nLayers)+1;
    tmatches = cell(nScales,1);
    tconfidence = zeros(nScales,1);
    parfor k=1:nScales
        des_1 = descriptors_1{Octave1(k),Layer1(k)};
        des_2 = descriptors_2{Octave2(k),Layer2(k)};
        [tmatches{k},tconfidence(k)] = Match_Keypoint(des_1,des_2,Error,K);
    end
    for k=1:nScales
        matches{Octave1(k),Layer1(k),Octave2(k),Layer2(k)} = tmatches{k};
        confidence(Octave1(k),Layer1(k),Octave2(k),Layer2(k)) = tconfidence(k);
    end
else
    if scl_flag
        for octave2=1:nOctaves2
            for octave1=1:nOctaves1
                for layer2=1:nLayers
                    for layer1=1:nLayers
        [matches{octave1,layer1,octave2,layer2},...
            confidence(octave1,layer1,octave2,layer2)] = Match_Keypoint(...
            descriptors_1{octave1,layer1},descriptors_2{octave2,layer2},Error,K);
                    end
                end
            end
        end
    else
        for octave=1:min(nOctaves1,nOctaves2)
            for layer2=1:nLayers
                for layer1=1:nLayers
        [matches{octave,layer1,octave,layer2},...
            confidence(octave,layer1,octave,layer2)] = Match_Keypoint(...
            descriptors_1{octave,layer1},descriptors_2{octave,layer2},Error,K);
                end
            end
        end
    end
end

%% Optimizing
Matches = cell(nOctaves1,nOctaves2);
Confidence = zeros(nOctaves1,nOctaves2);
for octave1=1:nOctaves1
    for octave2=1:nOctaves2
        matches_t = [];
        for layer1=1:nLayers
            for layer2=1:nLayers
                matches_t = [matches_t; matches{octave1,layer1,octave2,layer2}];
            end
        end
        if size(matches_t,1)>20
            matches_t = matches_t(:,[3:4,1:2,5,8:9,6:7,10]);  % Switch kps and kps_t
            [~,index1,~] = unique(matches_t(:,1:2),'rows');
            matches_t = matches_t(index1,:);
            [~,index2,~] = unique(matches_t(:,6:7),'rows');
            matches_t = matches_t(index2,:);
        end
        if size(matches_t,1)>20
            Matches{octave1,octave2} = matches_t;
            Confidence(octave1,octave2) = size(matches_t,1);
        end
    end
end
[max_O1,max_O2] = find(Confidence==max(max(Confidence)));

MMatches = [];
for i = 1-min(max_O1,max_O2):min(nOctaves1-max_O1,nOctaves2-max_O2)
    matches_t = Matches{max_O1+i,max_O2+i};
    if size(matches_t,1)>3
        MMatches = [MMatches; matches_t];
    end
end
[~,index1,~] = unique(MMatches(:,1:2),'rows');
MMatches = MMatches(index1,:);
[~,index2,~] = unique(MMatches(:,6:7),'rows');
MMatches = MMatches(index2,:);

%% One last outlier removal
NCMs = zeros(K,1); indexPairs = cell(K,1);
for k = 1:K
    [~,~,indexPairs{k}] = Outlier_Removal(MMatches(:,1:5),MMatches(:,6:end),Error);
    NCMs(k) = sum(indexPairs{k});
end
[~,maxIdx] = max(NCMs);
indexPairs = indexPairs{maxIdx};
MMatches = MMatches(indexPairs,:);
cor1 = MMatches(:,1:5); cor2 = MMatches(:,6:end);