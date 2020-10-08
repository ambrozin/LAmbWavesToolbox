figure(1) 
plot(real(arrayRes)')



S = fft(real(arrayRes(51:60,:)),[],2);

S2 = S(:,41);

Rxx = S2*S2'; 

[V D] = eig(Rxx);

plot(diag(D))
% plot(angle(V(end,:)))

dx = array.dx*1e-3; 
k = -500:0.1:500; 
ek = exp(1j*k.*[-4.5:4.5]'.*dx);

for i = 1: length(k); 
 
        P(i) = 1 ./ abs(sum(ek(:,i)'*V(:,1:9))).^2;

end

figure(3131)
plot(k,P)