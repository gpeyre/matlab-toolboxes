function out = cwpt2_interface(in, direction, filter_type, trans_type, btree)
% function out = cwpt2_interface(in, direction, filter_type, trans_type, btree)
%
% Forward 2D Continuous Wavelet Transform.
% Inverse 2D Continuous Wavelet Transform.
% Shifted Inverse 2D Continuous Wavelet Transform.
%
% Forward 2D Continuous Wavelet Packet Transform.
% Inverse 2D Continuous Wavelet Packet Transform.
%
% This is a wrapper for cwpt2, cwpt2i, c3wpt, c3wpti.
% It allows to select arbitrary filters to be used for the transform.
% It allows classical (quad) or 3-way wavelet elementary transforms.
%
% direction =   'forward' | 'inverse'
% filter_type = '7-9' | 'spline'
% trans_type =  'quad' | 'tri'
% btree =       a) for wavelet packet transforms, the corresponding btree
%               b) for forward wavelet transform, the number of scales
%               c) for inverse wavelet transform, must bet set to 0
%               d) for shifted inverse wavelet transform, must be set to 1
%
% A forward transform given a M x N image provides a M x N x S 3-D matrix.
%
% ================
% For Wavelet Packet Transform, type 'help cwpt2' for more information about the
% order of subbands.
%
% ================
% For Wavelet Transform, S = ( 1 + 3 * nb_scales ) and the order of subbands is:
%
% HL, LH, HH (scale 1),
% HL, LH, HH (scale 2),
% ...
% LL, HL, LH, HH (scale nb_scales),
%
% where the first filtering direction is the first memory raster direction of data.
% For instance, LH means for Matlab images low-filtering along columns (first dimension)
% and high-filtering along rows (second dimension).
%
% ================
% For 3-way Wavelet Transform, S = ( 1 + 2 * nb_scales ) and the order of subbands is:
%
% H0, 0H (scale 1),
% H0, 0H (scale 2),
% ...
% LL, H0, 0H (scale nb_scales),
%
% ================
% A shifted inverse (classical or 3-way) wavelet transform considers a transform of
% scale N as a transform of scale (N+1)
%

if ~isequal(class(in),'double')
    warning('Input value should be of type ''double''.');
    in = double(in);
end;

switch lower(filter_type)
case '7-9'
    % analysis
    alf = [ 0.03782845550699547033,-0.02384946501937996663,-0.1106244044184234859, 0.3774028556126538536, 0.8526986790094033264, 0.3774028556126538536,-0.1106244044184233194,-0.02384946501937999092, 0.03782845550699546339];
    ahf = [ 0.06453888262893847649,-0.04068941760955852721,-0.4180922732222120963, 0.7884856164056645023,-0.4180922732222124294,-0.04068941760955835374, 0.06453888262893842098];
    % reconstruction
    rlf = [-0.03226944131446923825,-0.02034470880477926361, 0.2090461366111060482, 0.3942428082028322511, 0.2090461366111062147,-0.02034470880477917687,-0.03226944131446921049];
    % ahf / 2; * -1 ^ odd
    rhf = [ 0.01891422775349773516, 0.01192473250968998331,-0.05531220220921174296,-0.1887014278063269268, 0.4263493395047016632,-0.1887014278063269268,-0.0553122022092116597, 0.01192473250968999546, 0.01891422775349773169];
    % alf / 2; * -1 ^ even
case 'spline'
    h = MakeONFilter_copy('Battle',3);
    g = ((-1 * ones(1,41)).^(0:40)).*h;

    alf = h;
    ahf = g;
    rlf = h./2;
    rhf = g./2;
otherwise
    error('Unknown filter type');
end;

switch lower(direction)
case 'forward'

    if length(btree)==1                         % wavelet; btree = nb_scales
        btree = cwpt2_btree(btree, 2);          % 2 = wavelet
    else                                        % packets
        if isequal(trans_type,'tri')
            warning('3-way transform has not been tested for WP transform');
        end;
    end;

    switch lower(trans_type)
    case 'quad'
        out = cwpt2(in, btree, alf, floor(length(alf)/2), ahf, floor(length(ahf)/2));
    case 'tri'
        out = c3wpt(in, btree, alf, floor(length(alf)/2), ahf, floor(length(ahf)/2));
    otherwise
        error('Unknown transform type');
    end;

case 'inverse'

    if length(size(in))~=3
        error('The input argument must be a M x N x S 3-D matrix.');
    end;

    if length(btree)==1                         % wavelet; btree = scale_shift
        if btree~=0 && btree~=1
            error('btree parameter must be 0 or 1 for inverse wavelet transform')
        end;

        nb_scales = size(in,3) - 1;
        if isequal(trans_type,'quad')
            if mod(nb_scales,3)~=0
                error('Invalid data: the size along the dimension 3 must be ( 1 + 3 * nb_scales ).');
            end;
            nb_scales = btree + (nb_scales / 3);
        else
            if mod(nb_scales,2)~=0
                error('Invalid data: the size along the dimension 3 must be ( 1 + 2 * nb_scales ).');
            end;
            nb_scales = btree + (nb_scales / 2);
        end;

        btree = cwpt2_btree(nb_scales, 2);      % 2 = wavelet
    else                                        % packets
        if isequal(trans_type,'tri')
            warning('3-way transform has not been tested for WP transform');
        end;
    end;

    switch lower(trans_type)
    case 'quad'
        out = cwpt2i(in, btree, rlf, floor(length(rlf)/2), rhf, floor(length(rhf)/2));
    case 'tri'
        out = c3wpti(in, btree, alf, floor(length(alf)/2), ahf, floor(length(ahf)/2), rlf, floor(length(rlf)/2), rhf, floor(length(rhf)/2));
    otherwise
        error('Unknown transform type');
    end;

otherwise
    error('Incorrect transform direction');
end;

% if the output argument is a string, this is an error message
if isequal(class(out),'char')
    error(out);
end;

% 2D Continuous Wavelet Packet Transform package
% (c) 2002-2005 Let It Wave, all rights reserved
