function residuals = reprojectionResidualStr(ObsIdx, ObsVal, px, py, f, Mot, Str)

% Variable to optimize for is only structure. So the motion is passed is as
% a constant into the function. Motion should be of the form 3x2xK while passing
% into the function

nCam = size(ObsIdx,1);

% Assuming that this function would be used only when motion is given
% and optimizing for str only. Meaning, x passed in should be having
% 3xnPts elements
if nargin ~= 7
    disp('Str only BA requires 7 arguments! Sth wrong. Abort!');
end

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