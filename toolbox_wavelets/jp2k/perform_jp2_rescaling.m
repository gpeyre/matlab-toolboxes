function MW = perform_jp2_rescaling(MW,Jmin,dir)

n = size(MW,1);
Jmax = log2(n)-1;

for j=Jmax:-1:Jmin
    q_min = 1;
    if j==Jmin
        q_min = 0;
    end
    for q=q_min:3
        [selx,sely] = compute_quadrant_selection(j,q);
        MW(selx,sely) = 2^(-(Jmax-j+1)*dir) * MW(selx,sely);
    end
end