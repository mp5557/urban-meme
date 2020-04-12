

classdef CubicSpline
    %CUBICSPLINE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        % A^T A x = A^T b
        ata_
        atb_
        del_x_
        num_ctr_
        res_x_
        param_
        A
        AInv
        ind_mat_
    end
    
    methods
        function obj = CubicSpline(del_x, num_ctr)
            %CUBICSPLINE 构造此类的实例
            %   此处显示详细说明
            obj.del_x_ = del_x;
            obj.num_ctr_ = num_ctr;
            obj.ata_ = zeros(2 * num_ctr);
            obj.atb_ = zeros(2 * num_ctr, 1);

            obj.A = [1 0 0 0
                     0 1 0 0
                     1 1 1 1
                     0 1 2 3];
            obj.AInv = inv(obj.A);
        end
        
        function obj = Update(obj,x, y)
            %METHOD1 此处显示有关此方法的摘要
            %   此处显示详细说明
            for ii = 1 : length(x)
                [x_base, x_vec] = obj.GetBase(x(ii));
                a = zeros(2*obj.num_ctr_, 1);
                a(x_base*2+1:x_base*2+4) = x_vec * obj.AInv;
                obj.ata_ = obj.ata_ + a * a';
                obj.atb_ = obj.atb_ + y(ii) * a;
            end
        end
        
        function [x_base, x_vec] = GetBase(obj, x)
            x_norm = x / obj.del_x_;
            x_base = Clamp(floor(x_norm), 0, obj.num_ctr_ - 2);
            x_scale = x_norm - x_base;
            x_vec = x_scale .^ [0,1,2,3];
        end
        
        function obj = Solve(obj)
            obj.res_x_ = obj.ata_ \ obj.atb_;
            for ii = 0 : obj.num_ctr_ - 2
                obj.param_(:, ii+1) = obj.AInv *  obj.res_x_(2*ii+1:2*ii+4);
            end
        end
        
        function [y] = Eval(obj, x)
            y = zeros(size(x));
            for ii = 1 : length(x)
                [x_base, x_vec] = obj.GetBase(x(ii));
                y(ii) = dot(obj.param_(:, x_base+1), x_vec);
            end
        end
    end
end

function [clamped] = Clamp(x, vmin, vmax)
fmin = x < vmin;
fmax = x > vmax;
clamped = fmin * vmin + fmax * vmax + (~fmin && ~fmax) * x;
end