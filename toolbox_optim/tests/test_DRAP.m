clear all
clf

projV=@(x,A)A*inv(A'*A)*A'*x;
V1=[0 1]';
V2=[1 1]';

% Alternating projections.
x=[0.5 1]';
subplot(121)
plot([0 V1(1)],[0 V1(2)],'-b','LineWidth',2);hold on
plot([0 V2(1)],[0 V2(2)],'-b','LineWidth',2);
plot(x(1),x(2),'.k','MarkerSize',20);
for i=1:5
  xV1=projV(x,V1);
  h=arrow(x,xV1,'Length',5);
  set(h,'LineStyle',':','EdgeColor','r','FaceColor','r');
  text(xV1(1),xV1(2),['P_{V_1}(x)'])
  x=projV(xV1,V2);
  h=arrow(xV1,x,'Length',5);
  text(x(1),x(2),['P_{V_2}(P_{V_1}(x))'])
  plot(x(1),x(2),'.k','MarkerSize',20);
end
axis([0 1 0 1]);box on;axis equal;axis tight
title('Alternating projections');

% Average Alternating Reflections/DR.
x=[0.5 1]';
subplot(122)
plot([0 V1(1)],[-1 V1(2)],'-b','LineWidth',2);hold on
plot([-1 V2(1)],[-1 V2(2)],'-b','LineWidth',2);
plot(x(1),x(2),'.k','MarkerSize',20);
for i=1:5
  xV1=2*projV(x,V1)-x;
  h=arrow(x,xV1,'Length',5);
  set(h,'LineStyle',':','EdgeColor','r','FaceColor','r');
  text(xV1(1),xV1(2),['R_{V_1}(x)'])
  xV2=2*projV(xV1,V2)-xV1;
  h=arrow(xV1,xV2,'Length',5);
  set(h,'LineStyle',':','EdgeColor','r','FaceColor','r');
  text(xV2(1),xV2(2),['R_{V_2}(R_{V_1}(x))'])
  x=(xV2+x)/2;
  arrow(xV2,x,'Length',5);
  text(x(1),x(2),['(R_{V_2}(R_{V_1}(x))+x)/2'])
  plot(x(1),x(2),'.k','MarkerSize',20);
end
axis([-1 1 -1 1]);box on;axis equal;axis tight
title('Average Alternating Reflections/DR');
