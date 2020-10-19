function residuals = reprojectionResidualIntrinsics(ObsIdx,ObsVal,K,Mot,Str)

% (Constant) ObsIdx: index of KxN for N points observed by K cameras, sparse matrix
% (Constant) ObsVal: 2xM for M observations
% px,py: princple points in pixels
% f: focal length in pixels
% Mot: 3x2xK for K cameras
% Str: 3xN for N points

nCam = size(ObsIdx,1);

if nargin == 3
    [K, Mot, Str] = unpackKMotStr(nCam, K);
end

K   = reshape(K, 3, 3);
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
    % Apply the intrinsics
    homo_xy = K*[TRX;TRY;TRZ];
    % un-homogenize to get coordinates
    x = homo_xy(1,:)./homo_xy(3,:);
    y = homo_xy(2,:)./homo_xy(3,:);
    % Take the actual position of the point in the aime
    ox = ObsVal(1,validIdx);
    oy = ObsVal(2,validIdx);
    % Difference is the residual
    residuals = [residuals [x-ox; y-oy]];    
end

residuals = residuals(:);
