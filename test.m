z = -pi:0.1:pi;
y = (cos(z)-0.8)./((cos(z)).^2 + cos(z)+0.5);
plot(z, y);