function [xx, tt, uu] = hypermetesting(tspan, xspan, u0, ul, scheme, cfl, dx, dt)
  % funzione per risolvere le equazioni iperboliche di primo grado omogenee
  % tspan  - estremi intervallo temporale
  % xspan  - estremi intervallo spaziale
  % u0     - function handle condizione iniziale
  % ul     - function handle condizione al bordo
  % scheme - scelta del metodo esplicito
  % cfl    - parametro da rispettare per la stabilità
  % dx     - passo spaziale
  % dt     - passo temporale
  % 
  %
  lambda = dt/dx;
  a = cfl*lambda^(-1);
  xx = xspan(1):dx:xspan(2);
  tt = tspan(1):dt:tspan(2);
  uu = zeros(size(xx,2), size(tt,2));
  
  % condizione iniziale al tempo 0
  uu(:,1) = u0(xx);
  uu(1,:) = ul(tt);
  
  
  
  
%condizioni al contorno: valori della funzione u agli estremi ad ogni
%istante
%{
for k=1:length(xx)
    u(1,k) =ul(k);
end

%condizione iniziale
for j=1:length(tt)
    u(j,1) = u0(j);
end
%}
 %uu(1,1)=u0(1);
  
  if scheme == 1 % Lax Friedrichs   
    for n = 1:size(uu,2)-1
      for j = 2:size(uu,1)-1
            uu(j,n+1) = 0.5*(uu(j+1,n)+uu(j-1,n)) - 0.5*lambda*a*(uu(j+1,n)-uu(j-1,n));
      end
    end
    
  elseif scheme == 2 % Lax Wendroff 
    for n = 1:size(uu,2)-1
      for j = 2:size(uu,1)-1
        uu(j,n+1) = uu(j,n) - 0.5*lambda*a*(uu(j+1,n)-uu(j-1,n)) + 0.5*lambda^2*a^2*(uu(j+1,n)-2*uu(j,n) + uu(j-1,n));
      end
      end 
        
  elseif scheme == 3 % Upwind
    
    for n = 1:size(uu,2)-1
      for j = 2:size(uu,1)-1
        uu(j,n+1) = uu(j,n) - 0.5*lambda*a*(uu(j+1,n)-uu(j-1,n)) + 0.5*lambda*abs(a)*(uu(j+1,n)-2*uu(j,n)+uu(j-1,n));
      end
    end
   %elseif scheme == 4 % originale
    
   % for n = 1:size(uu,2)-1
   %   for j = 2:size(uu,1)-1
   %     uu(j,n) = uex_t(j*dx,n*dt);  
   %     if n==0
   %         uex_t(j,n)= sin(j); 
   %     end
   %     if j==0
   %         uex_t(j,n)= 0;
   %     end  
   %   end
   %   end
      end    
end