function [phi,out,Dice_history,Jaccard_history] =DRLSE_ADMM(Img,phi_0, g,gr, mu,alpha,beta, epsilon, r_u, r_q, r_varphi, r_p,iter) 
Dice_history=[];
Jaccard_history=[];


[Nx,Ny] = size(g);
% size(g)
% usufelu for FFT-based minimization
[Y,X] = meshgrid(0:Ny-1,0:Nx-1); 
auxFFT = cos(2*pi/Nx*X)+cos(2*pi/Ny*Y)-2;

% initialize Lagrange multipliers and splitting variables
lambda_u = zeros(Nx,Ny);
lambda_qx = zeros(Nx,Ny);
lambda_qy = zeros(Nx,Ny);
lambda_varphi = zeros(Nx,Ny);
lambda_px = zeros(Nx,Ny);
lambda_py = zeros(Nx,Ny);
phi=phi_0;
varphi = phi_0;
u = H(phi_0,epsilon);
q_x = Fx(u);
q_y = Fy(u);
 p_x = Fx(phi);
 p_y = Fy(phi);
% verb = 1; fig_it = figure;


for k=1:iter
%    
     % u^n+1 solved in Fourier domain(u-subproblem)
     b = r_u*(H(varphi,epsilon)- lambda_u/r_u) - beta*gr-r_q*( Bx(q_x+lambda_qx/r_q)+By(q_y+lambda_qy/r_q) );
      u = real(ifft2( fft2( b )./( -2*r_q*auxFFT + r_u )));
    
    % q^n+1 solved by shrinkage (q-subproblem)
    Fxu = Fx(u);
    Fyu = Fy(u);
    [q_x,q_y] = shrink(Fxu - lambda_qx/r_q,Fyu - lambda_qy/r_q,alpha*g/r_q);
    
     % varphi^n+1 approx initialization and Newton method
     % (varphi-subproblem)
    v = u + lambda_u/r_u;
    z = phi- lambda_varphi/r_varphi;
    E_0 = 0.5*r_varphi*z.^2 + 0.5*r_u*( 0.5 - v ).^2;
    E_z = 0.5*r_u*( H(z,epsilon) - v ).^2;
        % initialize with solution to ideal Heaviside
    cond_m2 = ( z>=0 ).*( v>=0.5 ) + ( z<0 ).*( v<0.5 );
    varphi = z.*cond_m2 + (1-cond_m2).*( E_z < E_0 ).*z;  
        % newton method to find zeros of derivative
    for j=1:3
         df = r_varphi*( varphi - z) + r_u*( H(varphi,epsilon) - v ).*D(varphi,epsilon);
         hf = r_varphi + r_u*D(varphi,epsilon).*D(varphi,epsilon) + r_u*(H(varphi,epsilon) - v ).*dD(varphi,epsilon);
         varphi = varphi - df./( hf + (abs(hf)<1e-3) ); 
    end
    
    % lambda^n+1 update lagrangians (update lagrange multiplers)
    lambda_u = lambda_u + r_u* (u-H(varphi,epsilon));
    lambda_qx =  lambda_qx + r_q* (q_x-Fxu);
    lambda_qy = lambda_qy + r_q* (q_y-Fyu);
    lambda_varphi = lambda_varphi + r_varphi* (varphi- phi);
    
    
        for jphi=1:2
              ex = Fx(phi)-lambda_px/r_p;
        ey = Fy(phi)-lambda_py/r_p;
        norm_e = sqrt( ex.^2 + ey.^2 );
        
         if and(ex,ey) ~=0
           p_x = (mu+r_p*norm_e./(mu+r_p)).*ex./(norm_e);
           p_y = (mu+r_p*norm_e./(mu+r_p)).*ey./(norm_e);
        elseif and(ex,ey) == 0
          p_x = mu+r_p/(mu+r_p);
          p_y = mu+r_p./(mu+r_p);
        else
            p_x=ex;
            p_y =ey;
        end
        end
    
    % phi^n+1 solved in Fourier domain (phi-subproblem)
        b = r_varphi* (varphi+lambda_varphi/r_varphi) -r_p.* ( Bx(p_x+lambda_px/r_p)+By(p_y+lambda_py/r_p) );
        phi = real(ifft2( fft2( b )./( -2*r_p.*auxFFT + r_varphi )));
    
    
    lambda_px = lambda_px + r_p.*(p_x-Fx(phi));
    lambda_py = lambda_py + r_p.*(p_y-Fy(phi));
   out.iter =iter;


end
 return; 


% shrinkage function
function [xs,ys] = shrink(x,y,lambda)

s = sqrt(x.*conj(x)+y.*conj(y));
ss = s-lambda;
ss = ss.*(ss>0);

s = s+(s<lambda);
ss = ss./s;

xs = ss.*x;
ys = ss.*y;

return;

% smooth approximation of Heaviside function
function v = H(phi,epsilon)
v = 0.5*( 1 + 2/pi*atan(phi/epsilon) );
return;

% smooth approximation of Dirac distribution
function v = D(phi,epsilon)
v = 1/pi* epsilon./ (epsilon^2 + phi.^2);
return;

% smooth approximation of derivative of the Dirac distribution
function v = dD(phi,epsilon)
v = -2/pi* epsilon*phi./ (epsilon^2 + phi.^2).^2;


% Forward derivative operator on x with boundary condition u(:,1,:)=u(:,Nx,:)
function Fxu = Fx(u)
[Ny,Nx] = size(u);
Fxu = circshift(u,[0 -1])-u;
Fxu(:,Nx) = zeros(Ny,1);

% Forward derivative operator on y with boundary condition u(1,:,:)=u(Ny,:,:)
function Fyu = Fy(u)
[Ny,Nx] = size(u);
Fyu = circshift(u,[-1 0])-u;
Fyu(Ny,:) = zeros(1,Nx);

% Backward derivative operator on x with boundary condition Bxu(:,1)=u(:,1)
function Bxu = Bx(u)
[Ny,Nx] = size(u);
Bxu = u - circshift(u,[0 1]);
Bxu(:,1) = u(:,1);
Bxu(:,Nx) = -u(:,Nx-1);

% Backward derivative operator on x with boundary condition Bxu(1,:)=u(1,:)
function Byu = By(u)
[Ny,Nx] = size(u);
Byu = u - circshift(u,[1 0]);
Byu(1,:) = u(1,:);
Byu(Ny,:) = -u(Ny-1,:);

