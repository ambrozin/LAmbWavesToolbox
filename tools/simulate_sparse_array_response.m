function Res = simulate_sparse_array_response(array,target,exct,f_vec,k_a, k_s,varargin)

if nargin == 6
    varargin{1} = 'a0s0'; 
end

Res = zeros(size(array.u,1), size(array.X,2), exct.n_sampl);
for iT = 1 : length(target.x)
    for j = 1: size(array.u,1)
        % source to target disance
        d1 = sqrt((array.X(array.u(j,:)>0)-target.x(iT)).^2 + (array.Y(array.u(j,:)>0)-target.y(iT)).^2);
        for i = 1:size(array.X,2)
            
            % target to detector disance
            d2 = sqrt((array.X(i)-target.x(iT)).^2+(array.Y(i)-target.y(iT)).^2);
            % total distance
            d = d1 + d2;
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
            
            %         resp_a = dispResponse(f_vec,k_a,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
            %         %         resp_s = dispResponse(f_vec,k_s,exct.fs,exct.n_sampl,exct.f_exc,exct.n_cycl,d);
            Res(j,i,:) =   squeeze(Res(j,i,:)) + target.reflectivity(iT)*hilbert(resp_a )';
            clc, disp(['Simulating array response: ' num2str(i/size(array.X,2)*100), '%'])
        end
    end
end