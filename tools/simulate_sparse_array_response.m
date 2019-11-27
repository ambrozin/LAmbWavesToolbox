function Res = simulate_sparse_array_response(array,target,exct,f_vec,k_a, k_s)

Res = zeros(size(array.u,1), size(array.X,2), exct.n_sampl);

for j = 1: size(array.u,1)
    % source to target disance
    d1 = sqrt((array.X(array.u(j,:)>0)-target.x).^2 + (array.Y(array.u(j,:)>0)-target.x).^2);
    for i = 1:size(array.X,2)
        
        % target to detector disance
        d2 = sqrt((array.X(i)-target.x).^2+(array.Y(i)-target.y).^2);
        % total distance 
        d = d1 + d2;
        resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        %         resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        Res(j,i,:) =  hilbert(resp_a );
        clc, disp(['Simulating array response: ' num2str(i/size(array.X,2)*100), '%'])
    end
end