function [K, Mot, Str] = unpackKMotStr(nCam, vec)

motion_cut = 3*2*nCam+5; % Mot variables at cut1+cut2

% K terms are fx, fy, cx, cy, s in that order in vec
% 
% fx = vec(1);
% fy = vec(2);
% cx = vec(3);
% cy = vec(4);
% s = vec(5);

K = [vec(1),vec(5),vec(3);0,vec(2),vec(4);0,0,1];
Mot = reshape(vec(6: motion_cut),3,2,[]);
Str = reshape(vec(motion_cut+1:end),3,[]);