[h, T_aw, i] = Nozzle_Flow();
T_i = 298;
c = 0.715; % J/(g*K)
rho = 1.81*100^3; % g/m^3
k = 92.67; % W/(m*K)
alpha = k/(rho*c); % m^2/s

time = [1,5,10,15,20,25,30,35,40,45,50,55,60];
position = [0,.02,.04,.06,.08,.10];
h_t = h(i);
T_aw_t = T_aw(i);
[T,X] = meshgrid(time,position);

ratio = @(T,X) 1-erf(X./(2*(alpha*T).^.5))-exp(h_t*X/k+h_t^2*alpha*T/k^2).*(1-erf(X./(2*(alpha*T).^.5)+h_t*(alpha*T).^.5/k));
Temperature = T_i + (T_aw_t-T_i)*ratio(T,X);

surf(T,X,Temperature)
xlabel("time (s)")
ylabel("radial depth (m)")
zlabel("Temperature (K)")