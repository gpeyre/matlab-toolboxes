function [Mi,err] = perform_linear_inpainting(M,type,mask,options)

% perform_linear_inpainting - perform inpainting with linear interpolation
%
%   [Mi,err] = perform_linear_inpainting(M,type,mask,options);
%
%   Depending on type, the solution is computed as
%       type='supperresol' : with bilinear interpolation
%       type='checkboard' and options.qsub==2 : quincunx interpolation
%       any other : heat diffusion, that minimise \int |grad(f)|^2.
%
%   Copyright (c) 2008 Gabriel Peyre

options.null = 0;
M0 = getoptions(options, 'M0', []);

err = [];
Mi = []; 
if strcmp(type, 'superresol')
    % linear interpolation
    qsub = getoptions(options, 'qsub', 1,1);
    [X,Y] = meshgrid(0:qsub:n-1+qsub,0:qsub:n-1+qsub);    
    [Xi,Yi] = meshgrid(0:n-1,0:n-1);
    Mz = M(1:qsub:end,1:qsub:end);
    Mz(:,end+1) = Mz(:,end); Mz(end+1,:) = Mz(end,:);
    Mi = interp2(X,Y,Mz,Xi,Yi);
end
if strcmp(type, 'checkboard') && qsub==2
    vm = getoptions(options, 'quincunx_vm', 6);
    Mi = perform_quincunx_interpolation(M,vm);
end
if isempty(Mi)
    % perform diffusion to compute initial starting point
    Mi = M;
    Mi(mask==1) = mean(M(mask==0));
    niterlin = getoptions(options, 'niterlin', 100);
    Slist = linspace(8,0,niterlin);
    options.bound = 'per';
    for i=1:niterlin
        progressbar(i,niterlin);
        if Slist(i)<4
            options.bound = 'sym';
        end
        Mi = perform_blurring(Mi,Slist(i),options);
        Mi(mask==0) = M(mask==0);   
        if not(isempty(M0))
            err(i) = snr(M0,Mi);
        end
    end
end