

nifti = niftiread('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/nifti/rpbm/1_AX_AdvDiff_Delta0030ms_b0_500_AX_AdvDiff_Delta0030ms_b0_500.nii.gz');
bvec=dlmread('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/nifti/rpbm/1_AX_AdvDiff_Delta0030ms_b0_500_AX_AdvDiff_Delta0030ms_b0_500.bval');
bval=dlmread('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/nifti/rpbm/1_AX_AdvDiff_Delta0030ms_b0_500_AX_AdvDiff_Delta0030ms_b0_500.bvec');
bmat = load('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/bmatrix/0030ms.txt');

bvec2=dlmread('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/nifti/rpbm/1_AX_AdvDiff_Delta0030ms_b0_500_AX_AdvDiff_Delta0030ms_b0_500_corr.bval');
bval2=dlmread('/Users/gabriellebaxter/Documents/AGING/data/S01/MRI/nifti/rpbm/1_AX_AdvDiff_Delta0030ms_b0_500_AX_AdvDiff_Delta0030ms_b0_500_corr.bvec');

% 
% nifti = niftiread('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016//derivatives/all/0050ms/dwiec.nii');
% bvec=dlmread('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016//derivatives/all/0050ms/dwiec.bvec');
% bval=dlmread('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016//derivatives/all/0050ms/dwiec.bval');
% bmat = load('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016/bmatrix/0050ms.txt');
% 
% bvec2=dlmread('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016/bvec_bval/0050ms.bval');
% bval2=dlmread('/Users/gabriellebaxter/Documents/HEAL_Volunteers/WCMyofascial1016/bvec_bval/0050ms.bvec');


B  =[bvec;bval].'; 

jellecode(B,bmat)
% b7 = bvalvec_to_B7(B);
% b7_bvec_bmat = bvalvec_to_B7(bmat);
% B2  =[bvec2;bval2].'; 
% b72 = bvalvec_to_B7(B2);
% 
% grad_table = bmatrix_to_grad(bmat);
% 
% % A = grad_table(:,2:4);
% % B = B(:,2:4);
% 
% R = A \ B;   % least-squares solution
% B_pred = A * R;
% max_err = max(abs(B_pred(:) - B(:)));

function jellecode(grad,bmatrix)

    cnt=[1,2,2,1,2,1];
    ind = [ 1     1
        1     2
        1     3
        2     2
        2     3
        3     3];
    
    normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
    grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
    
    normgrad     = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
    grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
    s = size(grad);
    b = [ones([s(1) 1],class(grad)) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
    disp(b)
    
    b = [ ones([s(1) 1],class(bmatrix)) -bmatrix ];
    b(:, [3 4 6]) = 2*b(:, [3 4 6]);
    disp(b)
 

end

function B7 = bvalvec_to_B7(grad)
    s = size(grad);
    if s(2) == 4
        cnt=[1,2,2,1,2,1];
        ind = [ 1     1
            1     2
            1     3
            2     2
            2     3
            3     3];
        
        normgrad = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
        grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
        
        normgrad     = sqrt(sum(grad(:, 1:3).^2, 2)); normgrad(normgrad == 0) = 1;
        grad(:, 1:3) = grad(:, 1:3)./repmat(normgrad, [1 3]);
        s = size(grad);
        b = [ones([s(1) 1],class(grad)) -(grad(:,ones(1,6)*4).*grad(:,ind(1:6,1)).*grad(:,ind(1:6,2)))*diag(cnt)];
        disp(b)
    
    elseif s(2) == 6 % b-matrix
        b = [ ones([s(1) 1],class(grad)) -grad ];
        b(:, [3 4 6]) = 2*b(:, [3 4 6]);
        disp(b)
    end
    B7 = b;
end

function grad_table = bmatrix_to_grad(Bmats)
    % Assume Bmats is Nx6: [Bxx, Bxy, Bxz, Byy, Byz, Bzz]
    N = size(Bmats, 1);
    grad_table = zeros(N,4);  % Output Nx4 [gx, gy, gz, b]
    
    for i = 1:N
        % Reconstruct full 3x3 B-matrix
        B = [ Bmats(i,1), Bmats(i,4), Bmats(i,5);  % [Bxx, Bxy, Bxz]
          Bmats(i,4), Bmats(i,2), Bmats(i,6);  % [Bxy, Byy, Byz]
          Bmats(i,5), Bmats(i,6), Bmats(i,3)];
          
        % Eigen-decomposition to find principal direction
        [V,D] = eig(B);
        
        % Find largest eigenvalue (b-value)
        [bval, idx] = max(diag(D));
        
        % Corresponding eigenvector = gradient direction
        gvec = V(:,idx);
        
        % Ensure unit length
        gvec = gvec / norm(gvec);
        
        % Store in gradient table
        grad_table(i,:) = [bval, gvec'];

    end
end