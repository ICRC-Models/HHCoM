x = 0 : 0.1 : 10;
y = 0 : 0.1 : 10;
[X , Y] = meshgrid(x , y);
z = (sin(X) - cos(Y)) .* sin(0.3 * X);
[r , c] = find(z >= (max(unique(max(z))) .* 0.98));
figure()
mesh(X , Y , z)
% hold on
% plot3(x(c) , y(r) , z(r , c) , 'r*');