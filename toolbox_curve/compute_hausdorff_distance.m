function [d_rms,d_inf] = compute_hausdorff_distance(c1,c2, p)

% compute_hausdorff_distance - compute the symmetric RMS hausdorff distance btw 2 curves
%
%   [d_rms,d_inf] = compute_hausdorff_distance(c1,c2, p);
%
%   p is the number of sampling points (default: 100).
%   Set p<0 if you do not want resampling.
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    p = 300;
end
if size(c1,2)<size(c1,1)
    c1 = c1';
end
if size(c2,2)<size(c2,1)
    c2 = c2';
end

d_rms = 0;

% re-sample the curve using arc-length
if p>0
    c1 = perform_curve_resampling(c1, p, 'nbpts');
    c2 = perform_curve_resampling(c2, p, 'nbpts');    
end

d_rms = 0;
d_inf = 0;
for tmp=1:2
    d1_rms = 0;
    d1_inf = 0;
    for i=1:length(c1)
        x = c1(:,i);
        % find closest point to x
        D = compute_distance_to_points(c2,x);
        d1_rms = d1_rms + min(D);
        d1_inf = max(d1_inf, min(D));
    end
    d_rms = d_rms + sqrt(d1_rms/length(c1));
    d_inf = d_inf + sqrt(d1_inf);
    % swap c1 and c2 to symmetrize error
    c = c1;
    c1 = c2;
    c2 = c;
end
d_rms = d_rms/2;
d_inf = d_inf/2;


