function residuals = reprojectionResidual(ObsIdx,ObsVal,px,py,f,Mot,Str)

% (Constant) ObsIdx: index of KxN for N points observed by K cameras, sparse matrix
% (Constant) ObsVal: 2xM for M observations
% px,py: princple points in pixels
% f: focal length in pixels
% Mot: 3x2xK for K cameras
% Str: 3xN for N points


nCam = size(ObsIdx,1);

if nargin==5
    % When motion and structure and f are all being optimized, x contains
    % f, camera motion and structure which needs to be reshaped.
    [Mot,Str,f] = unpackMotStrf(nCam,f);
elseif nargin==6
    % When f is given, the structure and motion needs to be recovered.
    [Mot,Str]   = unpackMotStrf(nCam,Mot);
end

Mot = reshape(Mot,3,2,[]);
Str = reshape(Str,3,[]);

residuals = [];
for c=1:nCam
    validPts = ObsIdx(c,:)~=0;
    validIdx = ObsIdx(c,validPts);
    
    % Do the extrinsics first (Rotation using axis angle)
    RP = AngleAxisRotatePts(Mot(:,1,c), Str(:,validPts));
    % Translate after the rotation
    TRX = RP(1,:) + Mot(1,2,c);
    TRY = RP(2,:) + Mot(2,2,c);
    TRZ = RP(3,:) + Mot(3,2,c);
    % Un-homogenize the coordinates
    TRXoZ = TRX./TRZ;
    TRYoZ = TRY./TRZ;
    % Do the intrinsics on the un-homogenized coordinates to get projected
    % point from current R,t estimate
    x = f*TRXoZ + px;
    y = f*TRYoZ + py;
    % Take the actual position of the point in the aime
    ox = ObsVal(1,validIdx);
    oy = ObsVal(2,validIdx);
    % Difference is the residual
    residuals = [residuals [x-ox; y-oy]];    
end

residuals = residuals(:);
