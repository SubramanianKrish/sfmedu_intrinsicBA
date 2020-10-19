function err = printReprojectionError(graph)

% Process motion to get a 3x2xK representation
nCam=length(graph.frames);
Mot = zeros(3,2,nCam);
for camera=1:nCam
    Mot(:,1,camera) = RotationMatrix2AngleAxis(graph.Mot(:,1:3,camera));
    Mot(:,2,camera) = graph.Mot(:,4,camera);
end

Str = graph.Str;
f  = graph.f;

% Assuming principal point is at the center of image
% May need to change later!
px = 0;
py = 0;

if isfield(graph, 'K')
    residuals = reprojectionResidualIntrinsics(graph.ObsIdx,graph.ObsVal,graph.K,Mot,Str);
else
    residuals = reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
end

err = 2*sqrt(sum(residuals.^2)/length(residuals));

fprintf('current reprojection error = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));