function Res = simulate_array_response(X,Y,target,exct,f_vec,k_a, k_s)

Res = zeros(size(X,2), exct.n_sampl);
for i = 1:size(X,2)
        d = sqrt((X(i)-target.x).^2+(Y(i)-target.y).^2); 
        resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
%         resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        Res(i,:) =  hilbert(resp_a );
    clc, disp(['Simulating array response: ' num2str(i/size(X,2)*100), '%'])
end