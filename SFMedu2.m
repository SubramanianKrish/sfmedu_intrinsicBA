clc;
disp('SFMedu: Structrue From Motion for Education Purpose');
disp('Version 2 @ 2014');
disp('Written by Jianxiong Xiao (MIT License)');


%% set up things
clear;
close all;
addpath(genpath('matchSIFT'));
addpath(genpath('denseMatch'));
addpath(genpath('RtToolbox'));

visualize = false;

%% data

% frames.images{1}='temple/temple0272.png';
% frames.images{2}='temple/temple0275.png';
% frames.images{3}='temple/temple0278.png';
% frames.images{4}='temple/temple0281.png';
% frames.images{5}='temple/temple0282.png';

frames.images{1}='images/B21.jpg';
frames.images{2}='images/B22.jpg';
frames.images{3}='images/B23.jpg';
frames.images{4}='images/B24.jpg';
frames.images{5}='images/B25.jpg';

% frames.images{1}='husky/husky1.jpg';
% frames.images{2}='husky/husky2.jpg';
% frames.images{3}='husky/husky3.jpg';
% frames.images{4}='husky/husky4.jpg';
% frames.images{5}='husky/husky5.jpg';


%{
frames.images{1}='images/kermit000.jpg';
frames.images{2}='images/kermit001.jpg';
frames.images{3}='images/kermit002.jpg';
%}

%{
frames.images{1} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.15.54.jpg';
frames.images{2} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.15.59.jpg';
frames.images{3} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.02.jpg';
frames.images{4} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.05.jpg';
frames.images{5} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.09.jpg';
frames.images{6} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.16.jpg';
frames.images{7} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.20.jpg';
frames.images{8} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.23.jpg';
frames.images{9} ='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.26.jpg';
frames.images{10}='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.29.jpg';
frames.images{11}='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.33.jpg';
frames.images{12}='/Users/xj/Dropbox/Camera Uploads/2014-09-03 14.16.36.jpg';
%}


%% data loading
frames.length = length(frames.images);

try
    frames.focal_length = extractFocalFromEXIF(frames.images{1});
catch
end

if ~isfield(frames,'focal_length') || isempty(frames.focal_length)
    fprintf('Warning: cannot find the focal length from the EXIF\n');
    frames.focal_length = 719.5459; % for testing with the B??.jpg sequences
    
%     frames.focal_length = 1520.0; % temple images
    
    % husky focal length
    % iphone SE specs https://www.devicespecifications.com/en/model/7a423ad7
%     ccd_width_mm = 4.8;
%     focal_mm = 4.15;
%     res_x = 3024; % width pixels
%     frames.focal_length = res_x * (focal_mm / ccd_width_mm);
end


maxSize = 1024;
frames.imsize = size(imread(frames.images{1}));
if max(frames.imsize)>maxSize
    scale = maxSize/max(frames.imsize);
    frames.focal_length = frames.focal_length * scale;
    frames.imsize = size(imresize(imread(frames.images{1}),scale));
end


frames.K = f2K(frames.focal_length);
disp('intrinsics:');
disp(frames.K);

%% SIFT matching and Fundamental Matrix Estimation
for frame=1:frames.length-1    
    % need to set this random seed to produce exact same result
    s = RandStream('mcg16807','Seed',10); 
    RandStream.setGlobalStream(s);
    
    % keypoint matching
    %pair = match2viewSURF(frames, frame, frame+1);
    pair = match2viewSIFT(frames, frame, frame+1);
    
    if visualize, showMatches(pair,frames); title('raw feature matching'); end
    
    if true % choose between different ways of getting E
        % Estimate Fundamental matrix
        pair = estimateF(pair);    
        % Convert Fundamental Matrix to Essential Matrix
        pair.E = frames.K' * pair.F * frames.K; % MVG Page 257 Equation 9.12
    else
        % Estimate Essential Matrix directly using 5-point algorithm
        pair = estimateE(pair,frames); 
    end
    

    if visualize, showMatches(pair,frames); title('inliers'); end

    % Get Poses from Essential Matrix
    pair.Rt = RtFromE(pair,frames);
    
    % Convert the pair into the BA format
    Graph{frame} = pair2graph(pair,frames);
    
    % re-triangulation
    Graph{frame} = triangulate(Graph{frame},frames);
    if visualize, visualizeGraph(Graph{frame},frames); title('triangulation'); end
    
    % outlier rejection
    % Graph{frame} = removeOutlierPts(Graph{frame});
    
    adjustFocalLength = false;
    adjustIntrinsics = false;
    
    % bundle adjustment for str, motion [original function call]
    Graph{frame} = bundleAdjustment(Graph{frame}, adjustFocalLength, adjustIntrinsics);
    
    % BA for str, motion and intrinsics
    % if a K of 3x3 is not present in the graph, BA will make one and
    % populate it
%     Graph{frame} = bundleAdjustment(Graph{frame}, adjustFocalLength, adjustIntrinsics);
    
    if visualize, visualizeGraph(Graph{frame},frames); title('after two-view bundle adjustment'); end
end


%% merge the graphs - standard way of merging
%close all
fprintf('\n\nmerging graphs....\n');

mergedGraph = Graph{1};
graph_error = [];
for frame=2:frames.length-1
    % merge graph
    
    % if an 3x3xnCams K is present in the two graphs, then the first graphs
    % K are taken for all the A_unique frames and AB_common frames. Meaning
    % the K estimate of a common frame in graph B is ignored. This is done
    % as we are anyways doing another BA after the merging at each stage of
    % merging. This therefore would only change the initial seed of the
    % optimization for the intrinsics to those belonging to A or the first
    % time that camera is seen in the merge process.
    
    mergedGraph = merge2graphs(mergedGraph,Graph{frame});
    
    % re-triangulation
    mergedGraph = triangulate(mergedGraph,frames);
    if visualize, visualizeGraph(mergedGraph,frames); title('triangulation'); end
    
    % outlier rejection
    % mergedGraph = removeOutlierPts(mergedGraph,10);
    
    % bundle adjustment
    mergedGraph = bundleAdjustment(mergedGraph, adjustFocalLength, adjustIntrinsics);
    
    % outlier rejection
    mergedGraph = removeOutlierPts(mergedGraph, 10);
    
    % bundle adjustment
    mergedGraph = bundleAdjustment(mergedGraph, adjustFocalLength, adjustIntrinsics);
    
    graph_error = [graph_error printReprojectionError(mergedGraph)];
    
    if visualize, visualizeGraph(mergedGraph,frames); title('after bundle adjustment'); end
end
% plot(graph_error);
%% Smart merging based on maximum common features

% fprintf('\n\nmerging graphs SMARTLY....\n');
% 
% % The starting point is the graph with most features
% max_matches = -1;
% merge_start_index = -1;
% for i=1:frames.length-1
%     if size(Graph{i}.matches,2) > max_matches
%         max_matches = size(Graph{i}.matches,2);
%         merge_start_index = i;
%     end
% end
% 
% mergedGraph = Graph{merge_start_index};
% current_index = merge_start_index;
% mergelog = 1:frames.length - 1;
% 
% while length(mergedGraph.frames) ~= frames.length
%     % ensures all graphs are merged. Code breaks if there exists graphs
%     % with no common frames
%     next_common = -1;
%     prev_common = -1;
%     if(current_index+1 <= length(mergelog))
%         next_graph_id = mergelog(current_index+1);
%         commonFrames_next = intersect(mergedGraph.frames,Graph{next_graph_id}.frames);
%         % We're gonna see whether the graph before or after the current
%         % graph shares maximum feature points with the current graph. We
%         % only check the number of commmon features in the first common
%         % frame found between the two graphs.
%         cameraIDA = find(mergedGraph.frames==commonFrames_next(1));
%         cameraIDB = find(Graph{next_graph_id}.frames==commonFrames_next(1));
%         
%         trA = find(mergedGraph.ObsIdx(cameraIDA,:)~=0);
%         xyA = mergedGraph.ObsVal(:,mergedGraph.ObsIdx(cameraIDA,trA));
% 
%         trB = find(Graph{next_graph_id}.ObsIdx(cameraIDB,:)~=0);
%         xyB = Graph{next_graph_id}.ObsVal(:,Graph{next_graph_id}.ObsIdx(cameraIDB,trB));
%         [xyCommon,iA,iB] = intersect(xyA',xyB','rows');
%         next_common = size(xyCommon, 1);
%     end
%     
%     if(current_index-1 > 0)
%         prev_graph_id = mergelog(current_index-1);
%         commonFrames_prev = intersect(mergedGraph.frames,Graph{prev_graph_id}.frames);
%         % We're gonna see whether the graph before or after the current
%         % graph shares maximum feature points with the current graph. We
%         % only check the number of commmon features in the first common
%         % frame found between the two graphs.
%         cameraIDA = find(mergedGraph.frames==commonFrames_prev(1));
%         cameraIDB = find(Graph{prev_graph_id}.frames==commonFrames_prev(1));
%         
%         trA = find(mergedGraph.ObsIdx(cameraIDA,:)~=0);
%         xyA = mergedGraph.ObsVal(:,mergedGraph.ObsIdx(cameraIDA,trA));
% 
%         trB = find(Graph{prev_graph_id}.ObsIdx(cameraIDB,:)~=0);
%         xyB = Graph{prev_graph_id}.ObsVal(:,Graph{prev_graph_id}.ObsIdx(cameraIDB,trB));
%         [xyCommon,iA,iB] = intersect(xyA',xyB','rows');
%         prev_common = size(xyCommon, 1);
%     end
%     
%     if next_common > prev_common
%         to_merge_id = current_index+1;
%         mergedGraph = merge2graphs(mergedGraph,Graph{mergelog(to_merge_id)});
%     elseif next_common < prev_common
%         to_merge_id = current_index-1;
%         mergedGraph = merge2graphs(mergedGraph,Graph{mergelog(to_merge_id)});
%     else
%         disp("reached end! Both next and prev common are -1");
%     end
%     
%     % update mergelog based on to_merge_id
%     mergelog([current_index, to_merge_id]) = [-1,-1];
%     h = 1;
%     templog = [];
%     while(h <= length(mergelog))
%         % squash those two -1s into one -1
%         if mergelog(h) == -1
%             h = h+1;
%             templog = [templog, -1];
%             current_index = length(templog);
%         else
%             templog = [templog, mergelog(h)];
%         end
%         h = h+1;
%     end
%     mergelog = templog;
%     
%     % re-triangulation
%     mergedGraph = triangulate(mergedGraph,frames);
%     if visualize, visualizeGraph(mergedGraph,frames); title('triangulation'); end
%     
%     % outlier rejection
%     % mergedGraph = removeOutlierPts(mergedGraph,10);
%     
%     % bundle adjustment
%     mergedGraph = bundleAdjustment(mergedGraph, adjustFocalLength, adjustIntrinsics);
%     
%     % outlier rejection
%     mergedGraph = removeOutlierPts(mergedGraph, 10);
%     
%     % bundle adjustment
%     mergedGraph = bundleAdjustment(mergedGraph, adjustFocalLength, adjustIntrinsics);    
%     
%     if visualize, visualizeGraph(mergedGraph,frames); title('after bundle adjustment'); end
% end
%%

%{
% outlier rejection
mergedGraph = removeOutlierPts(mergedGraph);

% bundle adjustment
mergedGraph = bundleAdjustment(mergedGraph);
if visualize, visualizeGraph(mergedGraph,frames); title('after bundle adjustment'); end

% bundle adjustment with focal length changes
mergedGraph = bundleAdjustment(mergedGraph,true);
if visualize, visualizeGraph(mergedGraph,frames); title('after bundle adjustment with focal length'); end
%}
%%
printReprojectionError(mergedGraph); % [for homework]

visualizeReprojection(mergedGraph,frames); % [for homework]
%%
points2ply('sparse.ply',mergedGraph.Str);

if frames.focal_length ~= mergedGraph.f
    disp('Focal length is adjusted by bundle adjustment');
    frames.focal_length = mergedGraph.f;
    frames.K = f2K(frames.focal_length);
    disp(frames.K);
end

if isfield(mergedGraph, 'K')
    % Populate frames with K if merged graph has a K created in BA
    frames.K = mergedGraph.K;
else
    % to not have another if statement in dense recon, we do same K
    % here
    frames.K = f2K(frames.focal_length);
end

%% dense matching

fprintf('dense matching ...\n');
for frame=1:frames.length-1
    Graph{frame} = denseMatch(Graph{frame}, frames, frame, frame+1);
end


%% dense reconstruction
fprintf('triangulating dense points ...\n');
for frame=1:frames.length-1
    clear X;
    P{1} = frames.K * mergedGraph.Mot(:,:,frame);
    P{2} = frames.K * mergedGraph.Mot(:,:,frame+1);
    %par
    for j=1:size(Graph{frame}.denseMatch,2)
        X(:,j) = vgg_X_from_xP_nonlin(reshape(Graph{frame}.denseMatch(1:4,j),2,2),P,repmat([frames.imsize(2);frames.imsize(1)],1,2));
    end
    X = X(1:3,:) ./ X([4 4 4],:);
    x1= P{1} * [X; ones(1,size(X,2))];
    x2= P{2} * [X; ones(1,size(X,2))];
    x1 = x1(1:2,:) ./ x1([3 3],:);
    x2 = x2(1:2,:) ./ x2([3 3],:);
    Graph{frame}.denseX = X;
    Graph{frame}.denseRepError = sum(([x1; x2] - Graph{frame}.denseMatch(1:4,:)).^2,1);
    
    Rt1 = mergedGraph.Mot(:, :, frame);
    Rt2 = mergedGraph.Mot(:, :, frame+1);
    C1 = - Rt1(1:3, 1:3)' * Rt1(:, 4);
    C2 = - Rt2(1:3, 1:3)' * Rt2(:, 4);
    view_dirs_1 = bsxfun(@minus, X, C1);
    view_dirs_2 = bsxfun(@minus, X, C2);
    view_dirs_1 = bsxfun(@times, view_dirs_1, 1 ./ sqrt(sum(view_dirs_1 .* view_dirs_1)));
    view_dirs_2 = bsxfun(@times, view_dirs_2, 1 ./ sqrt(sum(view_dirs_2 .* view_dirs_2)));
    Graph{frame}.cos_angles = sum(view_dirs_1 .* view_dirs_2);
    
    c_dir1 = Rt1(3, 1:3)';
    c_dir2 = Rt2(3, 1:3)';
    Graph{frame}.visible = (sum(bsxfun(@times, view_dirs_1, c_dir1)) > 0) & (sum(bsxfun(@times, view_dirs_2, c_dir2)) > 0);
end

% visualize the dense point cloud
if visualize
    figure
    for frame=1:frames.length-1
        hold on
        goodPoint =  Graph{frame}.denseRepError < 0.05;
        plot3(Graph{frame}.denseX(1,goodPoint),Graph{frame}.denseX(2,goodPoint),Graph{frame}.denseX(3,goodPoint),'.b','Markersize',1);
    end
    hold on
    plot3(mergedGraph.Str(1,:),mergedGraph.Str(2,:),mergedGraph.Str(3,:),'.r')
    axis equal
    title('dense cloud')
    for i=1:frames.length
        drawCamera(mergedGraph.Mot(:,:,i), frames.imsize(2), frames.imsize(1), frames.K(1,1), 0.001,i*2-1);
    end
    axis tight
end

% output as ply file to open in Meshlab (Open Software available at http://meshlab.sourceforge.net )
plyPoint = [];
plyColor = [];
for frame=1:frames.length-1
    goodPoint =  (Graph{frame}.denseRepError < 0.05) & (Graph{frame}.cos_angles < cos(5 / 180 * pi)) & Graph{frame}.visible;
    X = Graph{frame}.denseX(:,goodPoint);
    % get the color of the point
    P{1} = frames.K * mergedGraph.Mot(:,:,frame);
    x1= P{1} * [X; ones(1,size(X,2))];
    x1 = round(x1(1:2,:) ./ x1([3 3],:));
    x1(1,:) = frames.imsize(2)/2 - x1(1,:);
    x1(2,:) = frames.imsize(1)/2 - x1(2,:);
    indlin = sub2ind(frames.imsize(1:2),x1(2,:),x1(1,:));
    im = imresize(imread(frames.images{frame}),frames.imsize(1:2));
    imR = im(:,:,1);
    imG = im(:,:,2);
    imB = im(:,:,3);
    colorR = imR(indlin);
    colorG = imG(indlin);
    colorB = imB(indlin);
    plyPoint = [plyPoint X];
    plyColor = [plyColor [colorR; colorG; colorB]];
end

points2ply('dense.ply',plyPoint,plyColor);

fprintf('SFMedu is finished.\n Open the result dense.ply in Meshlab (Open Software available at http://meshlab.sourceforge.net ).\n Enjoy!\n');

