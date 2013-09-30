% Test for bezier curves

N = 1024;
cont = true;
hold on;
axis([0 1 0 1]);
C = [];
while cont
    [x,y,b] = ginput(1);
    cont = b==1;
    if cont
        plot(x,y, '.');
        axis([0 1 0 1]);
        C(:,end+1) = [x;y];
    end
end
compute_bezier_curve(C,N);
hold off;