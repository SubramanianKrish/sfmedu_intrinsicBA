function visualization = visualizeReprojection(graph, frames)

% Make local copies to process throught this function
% Done only as a backup to ensure original passed matrices aren't changed
% even by accident
ObsIdx = graph.ObsIdx;
ObsVal = graph.ObsVal;

% Get the number of camera views
nCam = size(ObsIdx,1);

% Assuming principal point is at the center of the image
% May need to change later! <TODO>!
px = 0;
py = 0;

f = graph.f;
Str = graph.Str;
% Convert motion to axis angle format to reuse existing code
Mot = zeros(3,2,nCam);
for camera=1:nCam
    Mot(:,1,camera) = RotationMatrix2AngleAxis(graph.Mot(:,1:3,camera));
    Mot(:,2,camera) = graph.Mot(:,4,camera);
end

for c=1:nCam
    validPts = ObsIdx(c,:)~=0;
    validIdx = ObsIdx(c,validPts);
    
    % Do the extrinsics first (Rotation using axis angle)
    RP = AngleAxisRotatePts(Mot(:,1,c), Str(:,validPts));
    
    % Translate after the rotation
    TRX = RP(1,:) + Mot(1,2,c);
    TRY = RP(2,:) + Mot(2,2,c);
    TRZ = RP(3,:) + Mot(3,2,c);
    if isfield(graph, 'K')
        % Apply the intrinsics
        homo_xy = graph.K*[TRX;TRY;TRZ];
        % un-homogenize to get coordinates
        x = homo_xy(1,:)./homo_xy(3,:);
        y = homo_xy(2,:)./homo_xy(3,:);
    
    else
        % Un-homogenize the coordinates
        TRXoZ = TRX./TRZ;
        TRYoZ = TRY./TRZ;

        % Do the intrinsics on the un-homogenized coordinates to get projected
        % point from current R,t estimate
        x = f*TRXoZ + px;
        y = f*TRYoZ + py;
    end 
    
    % Take the orignal position of the point in the image
    ox = ObsVal(1,validIdx);
    oy = ObsVal(2,validIdx);
    
    
    % Project non-observed points onto the iamge plane
    invalid_pts = (ObsIdx(c,:)==0);
    invalid_idx = ObsIdx(c, invalid_pts);
    rp_invalid = AngleAxisRotatePts(Mot(:,1,c), Str(:,invalid_pts));
    trx_invalid = rp_invalid(1,:) + Mot(1,2,c);
    try_invalid = rp_invalid(2,:) + Mot(2,2,c);
    trz_invalid = rp_invalid(3,:) + Mot(3,2,c);
    trxoz_invalid = trx_invalid./trz_invalid;
    tryoz_invalid = try_invalid./trz_invalid;
    x_invalid = f*trxoz_invalid + px;
    y_invalid = f*tryoz_invalid + py;
    
    % Plot the actual points
    image = imresize(imread(frames.images{c}),frames.imsize(1:2));
    figure
    imshow(image);
    hold on
    % Plot the projected x and y with a green + on the image
    plot(size(image,2)/2-x,size(image,1)/2-y,'g+');
    hold on
    % Plot the observed keypoint with a red x on the image
    plot(size(image,2)/2-ox,size(image,1)/2-oy,'rx');
    hold on
    % Plot the projected but unobserved keypoint with yello o in the image
     plot(size(image,2)/2-x_invalid,size(image,1)/2-y_invalid,'yo');
    hold on
    % plot the blue lines connecting the corresponding points
    x= [size(image,2)/2-x; size(image,2)/2-ox];
    y= [size(image,1)/2-y; size(image,1)/2-oy]; 
    plot(x,y,'-b');
    hold off 
end

visualization = true;
