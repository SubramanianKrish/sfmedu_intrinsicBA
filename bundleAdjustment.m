function graph = bundleAdjustment(graph, adjustFocalLength, adjustIntrinsics)

% convert from Rt matrix to AngleAxis

% Mot (refers to motion) is a 3x2x2 for our example
% first column is axis angle rep of rotation; second column is translation,
% each layer is one camera
nCam=length(graph.frames);
Mot = zeros(3,2,nCam);
for camera=1:nCam
    Mot(:,1,camera) = RotationMatrix2AngleAxis(graph.Mot(:,1:3,camera));
    Mot(:,2,camera) = graph.Mot(:,4,camera);
end

Str = graph.Str;
f  = graph.f;

% assume px, py=0 [answer to principal point question]
px = 0;
py = 0;

% bundle adjustment using lsqnonlin in Matlab (Levenberg-Marquardt)
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','off');

% Adjust structure [for homework]
% !!! fill in your code here
% Structure only BA
% residuals = reprojectionResidualStr(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
% fprintf('initial error before Structure only BA = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));
% [vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidualStr(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,x), Str(:),[],[],options);
% Str = reshape(vec,3,[]);
% fprintf('error after structure only BA = %f\n', 2*sqrt(resnorm/length(residuals)));

% Adjust motion [for homework]
% !!! fill in your code here

% Motion only BA
% residuals = reprojectionResidualMotion(graph.ObsIdx,graph.ObsVal,px,py,f,Str,Mot);
% fprintf('initial error before motion only BA = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));
% [vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidualMotion(graph.ObsIdx,graph.ObsVal,px,py,f,Str,x), Mot(:),[],[],options);
% Mot = reshape(vec,3,2,[]);
% fprintf('error after motion only BA = %f\n', 2*sqrt(resnorm/length(residuals)));

% adjust motion and structure
% Loss function for str AND mot BA

if ~adjustFocalLength && ~adjustIntrinsics
    residuals = reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
    fprintf('initial error = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));
    [vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,x), [Mot(:); Str(:)],[],[],options);
    [Mot,Str] = unpackMotStrf(nCam,vec);
    fprintf('error = %f\n', 2*sqrt(resnorm/length(residuals)));
end

if exist('adjustFocalLength','var') && adjustFocalLength
    % adjust focal length, motion and structure
    [vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,x), [f; Mot(:); Str(:)],[],[],options);
    [Mot,Str,f] = unpackMotStrf(nCam,vec);
    fprintf('error = %f\n', resnorm/length(residuals));
    graph.f = f;
end


% Optimize for full intrinsics(fx,fy,cx,cy,skew), motion and structure
if adjustIntrinsics && ~adjustFocalLength
    % Pass in f as a 2x1 matrix with [fx;fy]
    % px, py refer to principal points
    % stack K as a 3x3xnCams for future in graph
    if ~isfield(graph, 'K')
        graph.K = f2K(graph.f);
    end
    
    K = graph.K;
    
    residuals = reprojectionResidualIntrinsics(graph.ObsIdx, graph.ObsVal, K, Mot, Str);
    fprintf('initial error before IntrinicMotStr BA = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));
    K_optim_params = [K(1,1); K(2,2);K(1,3); K(2,3); K(1,2)];
    [vec,resnorm,residuals,exitflag] = lsqnonlin(@(x) reprojectionResidualIntrinsics(graph.ObsIdx,graph.ObsVal,x), [K_optim_params; Mot(:); Str(:)],[],[],options);
    [K, Mot,Str] = unpackKMotStr(nCam,vec);
    graph.K = K;
    fprintf('error after IntrinicMotStr BA = %f\n', 2*sqrt(resnorm/length(residuals)));
end

%residuals = reprojectionResidual(graph.ObsIdx,graph.ObsVal,px,py,f,Mot,Str);
%fprintf('final error = %f\n', 2*sqrt(sum(residuals.^2)/length(residuals)));

for camera=1:nCam
    graph.Mot(:,:,camera) = [AngleAxis2RotationMatrix(Mot(:,1,camera))  Mot(:,2,camera)];    
end
graph.Str = Str;
