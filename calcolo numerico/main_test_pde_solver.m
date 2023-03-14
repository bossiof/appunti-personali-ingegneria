  % --------- clear all -----------
  clear all
  close all
  clc

  % intervalli spaziali e temporali
  tspan = [0,1];
  xspan = [0,1];

   % altri parametri 
  a =0.1;                % velocità di propagazione
  cfl = 0.9;
  %[X,Y] = meshgrid(xx(1:end),tt);
  
  dt = 1.e-3;           % passo temporale
  dx = dt*a/cfl;        % passo spaziale
  
  xx = xspan(1):dx:xspan(2);
  tt = tspan(1):dt:tspan(2);
  
  % roba per plottare
  t = tt(5);            % istante temporale
  x_x = 200;
  
  % condizioni iniziali funzionanti
  %u0 = @(x)sin(x)+cos(0);
  %ul = @(t) sin(0)+cos(t);
  
  % condizioni iniziali modificate
  u0 = @(x) sin(pi*x);
  ul = @(t) 0*t;
 
  
  
    uex_t = @(x,t) (t+1).*sin(pi.*x);
  %if tt==0
  %   uex_t(xx,tt)= sin(xx); 
  %end
  %if xx==0
  %    uex_t(xx,tt)= 0;
  %end
  
  
  % -------- risoluzione dell'equazione con i 3 metodi espliciti stabili ------
  [xx, tt, uu_lf] = hypermetesting(tspan, xspan, u0, ul, 1, cfl, dx, dt );
  [xx, tt, uu_lw] = hypermetesting(tspan, xspan, u0, ul, 2, cfl, dx, dt );
  [xx, tt, uu_up] = hypermetesting(tspan, xspan, u0, ul, 3, cfl, dx, dt);
 %[xx, tt, uu_or] = hypermetesting(tspan, xspan, u0, ul, 4, cfl, dx, dt );

  
  
  %------- curve caratteristiche ---------
  
  fchar = @(t,x) a*t+x;
  x_char = linspace(-10,10,250);
  t_char = linspace(0,1,10);
  %figure(9)
  %for x = x_char
  %  hold on;
  %  plot(t_char, fchar(t_char,x))
  %endfor
  %xlim([0 0.5])
  %ylim([0 1])
  % title("Curve caratteristiche")
  %xlabel('x')
  %ylabel('t')

  
  
  
  % ------- soluzione esatta più soluzioni approssimate ---------
  
  figure(1)
  %{
  uex_t = @(xx,tt) tt.*sin(xx);
  if tt==0
     uex_t(xx,tt)= sin(xx); 
  end
  if xx==0
      uex_t(xx,tt)= 0;
  end
  %}
  hold on
  plot(xx, uex_t(xx,0.250), 'r');
  plot(xx, uu_lf(1:5001,0.250/dt),  'b--');
  plot(xx, uu_lw(1:5001,0.250/dt), 'm--');
  plot(xx, uu_up(1:5001,0.250/dt), 'g--');
  legend('esatta','lax freidrich','lax wendroff','upwind')
  hold off;
  
  % -------- errore nell'applicare i metodi ----------

  t_x = 1:size(tt,2);
  
  
  err_lf_re = abs(uex_t(xx(x_x), tt)-uu_lf(x_x,:))./abs(uex_t(xx(x_x), tt));
  err_lw_re = abs(uex_t(xx(x_x), tt)-uu_lw(x_x,:))./abs(uex_t(xx(x_x), tt));
  err_up_re = abs(uex_t(xx(x_x), tt)-uu_up(x_x,:))./abs(uex_t(xx(x_x), tt));
  
  err_lf = abs(uex_t(xx(x_x), tt)-uu_lf(x_x,:));
  err_lw = abs(uex_t(xx(x_x), tt)-uu_lw(x_x,:));
  err_up = abs(uex_t(xx(x_x), tt)-uu_up(x_x,:));
  
  figure(2)
  hold on
  plot(tt, err_lf,  'b--');
  plot(tt, err_lw, 'm--');
  plot(tt, err_up, 'g--');
  ylim([0 0.5])
  
  figure(3)
  subplot(3,1,1);
  hold on;
  plot(tt, uex_t(xx(x_x), tt), 'k');
  plot(tt, uu_lf(x_x,:),  'b--');
  plot(tt, uu_lw(x_x,:), 'm--');
  plot(tt, uu_up(x_x,:), 'g--');
  legend('esatta','lax freidrich','lax wendroff','upwind')
  hold off;
  
  subplot(3,1,2);
  hold on
  plot(tt, err_lf,  'b--');
  plot(tt, err_lw, 'm--');
  plot(tt, err_up, 'g--');
  legend('LF','LW','UP')
  title('errore')
  ylim([0 0.5])
  
  subplot(3,1,3);
  hold on;
  plot(tt, err_lf_re,  'b--');
  plot(tt, err_lw_re, 'm--');
  plot(tt, err_up_re, 'g--');
  legend('LF','LW','UP')
  title('errore relativo')
  
  
  
  %err_lf = norm(uex_t(xx(x_x), tt)-uu_lf(t_x),2);
  %err_lw = norm(uex_t(xx(x_x), tt)-uu_lw(t_x),2);
  %err_up = norm(uex_t(xx(x_x), tt)-uu_up(t_x),2);
  %figure(2)
  %hold on
  %plot(tt, err_lf,  'b--');
  %plot(tt, err_lw, 'm--');
  %plot(tt, err_up, 'g--');

  % vediamo se riusciamo a fare un grafico 3d
  
  % --------- grafici 3d -----------
  [X,Y] = meshgrid(xx(1:end),tt);
  figure(4)
  mesh(X,Y,uu_lw(1:end,:)')
  title('app. di lax wendroff con mesh')
   xlabel('asse x');
  ylabel('asse dei tempi t');
  zlabel('intensità u');
%{  
  figure(5)
  surf(X,Y,uu_lw(1:end,:)')
   title('app. di lax wendroff senza mesh')
%}  
  figure(6)
  mesh(X,Y,uu_lf(1:end,:)')
  title('app. di lax friedrichs con mesh')
   xlabel('asse x');
  ylabel('asse dei tempi t');
  zlabel('intensità u');
%{
  figure(7)
  surf(X,Y,uu_lf(1:end,:)')
  title('app. di lax friedrichs senza mesh')
 %}
  
  
  figure(8)
  mesh(X,Y,uu_up(1:end,:)')
  title('app. up-wind con mesh')
  xlabel('asse x');
  ylabel('asse dei tempi t');
  zlabel('intensità u');
  
 %{
 
  figure(9)
  surf(X,Y,uu_up(1:end,:)')
  title('app. up-wind senza mesh')
  
  %figure(6)
  %surf(X,Y,uex_t(X,Y))
  
  %figure(7)
  %mesh(X,Y,uex_t(X,Y))
  
  figure(10)
  surf(X,Y,abs(uex_t(X,Y)-uu_lw(1:end,:)'))
  title('errore assoluto rispetto a lax wendroff')
  
  figure(11)
  surf(X,Y,abs(uex_t(X,Y)-uu_lf(1:end,:)'))
  title('errore assoluto rispetto a lax friedrichs')
  
  figure(12)
  surf(X,Y,abs(uex_t(X,Y)-uu_up(1:end,:)'))
  title('errore assoluto rispetto ad upwind')
  
  figure(13)
  subplot(1,3,1)
  surf(X,Y,abs(uex_t(X,Y)-uu_lw(1:end,:)'))
  subplot(1,3,2)
  surf(X,Y,abs(uex_t(X,Y)-uu_lf(1:end,:)'))
  subplot(1,3,3)
  surf(X,Y,abs(uex_t(X,Y)-uu_up(1:end,:)'))
  title('insieme degli errori assoluti')
 %}
 
  
  figure(14);
  mesh(X,Y,uex_t(X,Y))
  xlabel('asse x');
  ylabel('asse dei tempi t');
  zlabel('intensità u');
  %surf(X,Y,uex_t(X,Y))
%{  
  figure(15)
  for i=1:1001 %// To see the results
    plot(xx,uu_lf(:,i))
    %plot(xx,uex_t(1:501,i))
    axis([-0.5 2.5 -3 3])
    drawnow
    pause(0.01)
end
%}
  %{
  figure(16);
  mesh(X,Y,uu_or(1:end,:)')
  title('figura originale campionata')
   xlabel('asse x');
  ylabel('asse dei tempi t');
  zlabel('intensità u');
  %}