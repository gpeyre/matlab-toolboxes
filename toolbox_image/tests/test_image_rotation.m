% test for image rotation

M = load_image('lena', 256);

theta = pi/50;
for i=1:100
    M = perform_image_rotation(M,theta);
    imageplot(M); drawnow;
end