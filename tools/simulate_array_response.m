function Res = simulate_array_response(X,Y,target,exct,f_vec,k_a, k_s, varargin)

Res = zeros(size(X,2), exct.n_sampl);
resp_a = zeros(1,exct.n_sampl); 
resp_s = zeros(1,exct.n_sampl); 


for j = 1 : length(target.x)
    for i = 1:size(X,2)
        d = sqrt((X(i)-target.x(j)).^2+(Y(i)-target.y(j)).^2);
       
        if nargin == 7 
            resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
            resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
        else 
            if strcmp(varargin{1}, 'a0')
                resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d); 
            end
            if strcmp(varargin{1}, 's0')
                resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d); 
            end
            if strcmp(varargin{1}, 'a0s0') % the same as default
                resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d); 
                resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d); 
            end
        end
        
        Res(i,:) =  Res(i,:) + target.reflectivity(j)*hilbert(resp_a + resp_s );
        clc, disp(['Simulating array response: ' num2str(i/size(X,2)*100), '%'])
    end
end