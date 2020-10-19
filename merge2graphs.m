function GraphAB = merge2graphs(GraphA,GraphB)


commonFrames = intersect(GraphA.frames,GraphB.frames);

[newFramesFromB,indexNewFramesFromB] = setdiff(GraphB.frames,GraphA.frames);

if isempty(commonFrames)
    GraphAB = [];
    return;
end

GraphAB = GraphA;

if isempty(newFramesFromB)
    return;
end


% add the non-overlapping frame first
firstCommonFrame = commonFrames(1);


% transform GraphB.Mot and GraphB.Str to be in the same world coordinate system of GraphA

% Each graph has it's own world coord system. Target is to get all of them
% under graph A's world coord system. Why first camera not Identity/world?
% Because BA done pairwise which moved the first camera around :3
% concatRts would essentially make two M matrices into one
% (just like T_AC = T_AB*T_BC transformation matrices)

RtBW2AW = concatenateRts(inverseRt(GraphA.Mot(:,:,GraphA.frames==firstCommonFrame)), GraphB.Mot(:,:,GraphB.frames==firstCommonFrame));
GraphB.Str = transformPtsByRt(GraphB.Str, RtBW2AW);
for i=1:length(GraphB.frames)
    GraphB.Mot(:,:,i) = concatenateRts(GraphB.Mot(:,:,i), inverseRt(RtBW2AW));
end

GraphAB.frames = [GraphA.frames newFramesFromB];
% Adding only the new frames changed M matrix as the common frames after
% the transformation to A world origin is same as the M of the common frame
% in A.
GraphAB.Mot(:,:,length(GraphA.frames)+1:length(GraphAB.frames)) = GraphB.Mot(:,:,indexNewFramesFromB);

% Add K just like motion if it exists
% if isfield(GraphA, 'K')
%     GraphAB.K(:,:,length(GraphA.frames)+1:length(GraphAB.frames)) = GraphB.K(:,:,indexNewFramesFromB);
% end
% No update on the K as we take the K from graph A.

% add the new tracks

for commonFrame = commonFrames
    
    cameraIDA = find(GraphA.frames==commonFrame);   cameraIDB = find(GraphB.frames==commonFrame);
    
    % trA refers to track A (which is set of indices(which can be used to
    % look up ground truth of the projection values) taken from a graph
    % where each line of obsIdx refers to one camera in graphA.
    
    % Obsval holds the original observed points in the 2xJ space where J is
    % the sum of all the original observed points from different images in that
    % graph. The indices got from trA can be used to access these values by
    % passing them into the column index. The rows here are just x and y.
    
    % Remember that ObsVal contains the coordinates of the corresponding
    % features between the 2 (or more) images stacked horizontally. hence
    % the 2xJ.
    
    % Althought the image is same between GraphA's common image and
    % graphB's common image, the feature points' corrdinates in their
    % ObsVal is different as the sift matching is done w.r.t another image.
    % For example im1, im2 in GA and im2, im3 in GB. ObsVal in GA contains
    % feature match coordinates between im1 and im2, while ObsVal in GB
    % contains feat. match coords betw im2 and im3. We seek the common
    % feature point found between both the image - xycommon
    
    trA = find(GraphA.ObsIdx(cameraIDA,:)~=0);
    xyA = GraphA.ObsVal(:,GraphA.ObsIdx(cameraIDA,trA));
    
    trB = find(GraphB.ObsIdx(cameraIDB,:)~=0);
    xyB = GraphB.ObsVal(:,GraphB.ObsIdx(cameraIDB,trB));
    
    % iA or iB is the index into xyA or xyB which results in that common
    % feature point.
    [xyCommon,iA,iB] = intersect(xyA',xyB','rows');
    xyCommon = xyCommon';
    
    % make the old track longer
    % Go through every common feature point
    for i=1:size(xyCommon,2)
        % Think of ObsIdx to be indexed as(cameraID, feature
        % feature point index)
        
        % Here idA refers to feature point index of that common feature in
        % grpahA. Feature point index means index into the ObsVal array.
        
        % iA holds the index of the common feature point (between graphA and graphB)
        % where A and B refer to the same image but different sets of 
        
        idA = trA(iA(i));
        idB = trB(iB(i));
        
        for j=1:length(indexNewFramesFromB)
            % Loop through every non-common frame in B
            % BObsIdx : look up the common feature point's index in the new
            % cameras of B. If non zero, means same feature seen in the new
            % camera view too (?)
            
            % find out the common feature point in the common frame. Look
            % up its correspondence in the non-common frames
            BObsIdx = GraphB.ObsIdx(indexNewFramesFromB(j),idB);
            
            % The BObsIdx would be zero if the correspondence for that
            % feature point did not exist in the non-common frame/view.
            
            % if correspondence exists in the non-common frame, go in
            if BObsIdx~=0
                % append the new feature point's location in the image  into merged
                % graph's feature point set
                GraphAB.ObsVal(:,end+1) = GraphB.ObsVal(:,BObsIdx);
                
                % Need to update the ObsIdx such that the cell indexed by
                % (camera ID, feature point index) needs to have the linear
                % index into the ObsVal matrix which would have the exact
                % pixel location (centralised coordinate system) of that
                % feature point.
                
                % ObsIdx helps look up n'th feature point coordinates in the 
                % i'th camera by holding "the linear index into the ObsVal
                % matrix" in a cell indexed by (cameraID (eg. i), local feature pointID(eg. n))
                
                % camera ID can be given by the number of already existing
                % cameras in A (this would've already included the common images/cameras)
                % plus the current camera in which we're
                % looking up the feature point. The local feature point
                % index ID is taken as the idA to maintain consistency for
                % later use. That way a column of merged graph's ObsIdx
                % would correspond to the same feature across all the
                % images :0
                GraphAB.ObsIdx(length(GraphA.frames)+j,idA) = size(GraphAB.ObsVal,2);            
            end
        end
    end
    
    % add the new tracks from common frame
    
    
    [xyNewFromB, iB] = setdiff(xyB',xyA','rows');
    xyNewFromB = xyNewFromB';
    
    for i=1:size(xyNewFromB,2)
        idB = trB(iB(i));
        
        GraphAB.ObsVal(:,end+1) = GraphB.ObsVal(:,GraphB.ObsIdx(cameraIDB,idB));
        GraphAB.ObsIdx(cameraIDA,end+1) = size(GraphAB.ObsVal,2);       
        GraphAB.Str(:,end+1) = GraphB.Str(:,idB);
        
        for j=1:length(indexNewFramesFromB)
            BObsIdx = GraphB.ObsIdx(indexNewFramesFromB(j),idB);
            if BObsIdx~=0
                GraphAB.ObsVal(:,end+1) = GraphB.ObsVal(:,BObsIdx);
                GraphAB.ObsIdx(length(GraphA.frames)+j,end) = size(GraphAB.ObsVal,2);            
            end
        end
    end    
    
end

% add the new tracks only among the completely new frames

newB = false(1,length(GraphB.frames));
newB(indexNewFramesFromB) = true;

tr2add = sum(GraphB.ObsIdx(~newB,:)~=0,1)==0 & sum(GraphB.ObsIdx(newB,:)~=0,1)>0;

if any(tr2add)
    
    ids = full(GraphB.ObsIdx(indexNewFramesFromB,tr2add));

    curValCnt = size(GraphAB.ObsVal,2);
    nonZerosID = find(ids(:)>0);

    GraphAB.ObsVal(:,(curValCnt+1):(curValCnt+length(nonZerosID))) = GraphB.ObsVal(:,ids(nonZerosID));

    idsNew = zeros(size(ids));
    idsNew(nonZerosID) = (curValCnt+1):(curValCnt+length(nonZerosID));

    GraphAB.ObsIdx(length(GraphA.frames)+1:end,size(GraphAB.ObsIdx,2)+1:size(GraphAB.ObsIdx,2)+size(idsNew,2)) = sparse(idsNew);

    GraphAB.Str(:,size(GraphAB.ObsIdx,2)+1:size(GraphAB.ObsIdx,2)+size(idsNew,2)) = GraphB.Str(:,tr2add);
end

return;



%{
    pointCount = size(MatchPairs{1}.matches,2);
    pointObservedValueCount = size(MatchPairs{1}.matches,2)*2;
    pointObservedValue(:,1:pointObservedValueCount) = [[MatchPairs{1}.matches(1:5,:) MatchPairs{1}.matches(6:10,:)]; -wTimePoints * ones(1,pointObservedValueCount)];
    pointObserved(1,1:pointCount)=1:pointCount;
    pointObserved(2,1:pointCount)=pointCount + (1:pointCount);
    previousIndex = 1:pointCount;
    
    pointCloud(:,1:pointCount) = MatchPairs{1}.matches(3:5,:);
    
    for frameID = 2:length(data.image)-1
        [~,iA,iB] = intersect(MatchPairs{frameID-1}.matches(6:7,:)',MatchPairs{frameID}.matches(1:2,:)','rows');
        
        
        alreadyExist = false(1,size(MatchPairs{frameID}.matches,2));
        alreadyExist(iB) = true;
        newCount = sum(~alreadyExist);
        
        
        currentIndex = zeros(1,size(MatchPairs{frameID}.matches,2));
        currentIndex(iB) = previousIndex(iA);
        currentIndex(~alreadyExist) = (pointCount+1):(pointCount+newCount);
        
        pointObservedValue(1:5,pointObservedValueCount+1:pointObservedValueCount+newCount+length(currentIndex)) = [MatchPairs{frameID}.matches(1:5,~alreadyExist) MatchPairs{frameID}.matches(6:10,:)];
        pointObservedValue(6,pointObservedValueCount+1:pointObservedValueCount+newCount+length(currentIndex)) = -wTimePoints;
        
        pointObserved(frameID  ,currentIndex(~alreadyExist)) = (pointObservedValueCount+1):(pointObservedValueCount+newCount);
        pointObservedValueCount = pointObservedValueCount + newCount;
        pointObserved(frameID+1,currentIndex) = (pointObservedValueCount+1):(pointObservedValueCount+length(currentIndex));
        pointObservedValueCount = pointObservedValueCount + length(currentIndex);
        
        
        pointCloud(:,pointCount+1:pointCount+newCount) = transformRT(MatchPairs{frameID}.matches(3:5,~alreadyExist), cameraRtC2W(:,:,frameID), false);
        
        pointCount = pointCount + newCount;
        
        previousIndex = currentIndex;
    end
%}