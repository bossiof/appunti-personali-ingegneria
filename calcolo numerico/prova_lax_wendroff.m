clear;

%parametri per definire l'advection equation e il range di spazio e tempo
Lmax = 1.0;
Tmax = 1.0;
c = 1;

 


%parametri per implementare il metodo di Lax-Wendroff
Nt = 1000; %numero di step temporali
Dt = Tmax/Nt; % singolo time step
Nx = 1000; %numero di step spaziali
Dx = Lmax/Nx; %singolo step spaziale
b = c*Dt/(2.*Dx); %parametro beta nell'implementazione alle differenze finite
x = 0:Dx:Lmax;
t = 0:Dt:Tmax;
[x,t] = meshgrid(x(1:end),t(1:end));
%il metodo di Lax Wendroff è stabile per b=<1/2 ma è diffuso a meno che
%abs(b)=1/2

%condizione iniziale: il valore iniziale della funzione u (ampiezza
%iniziale dell'onda)

%onda quadra:
uex_t = @(x,t) (t+1).*x.*2;
for i=1:(Nx+1)
            u(i,1)=uex_t(i,1);
        x(i)= (i-1)*Dx;
end


%condizioni al contorno: valori della funzione u agli estremi ad ogni
%istante
for i=1:(Nx+1)
    u(i,1)= sin(i/Nx);
    x(i)= (i-1)*Dx;
end
for j=1:(Nt+1)
    u(1,j)=0;
    t(j) = (j-1)*Dt; 
end
%implementazione del metodo di lax wendroff
for k=1:Nt  %ciclo temporale
    for i=2:Nx  %ciclo spaziale
        u(i,k+1) = 0.5*(u(i+1,k)+u(i-1,k))+b*(u(i+1,k)-u(i+1,k));
    end
end

mesh(x,t,u);
%axis([-0.5 1.5 -1 1.5 -1 2500])
%rappresentazione grafica dell'evoluzione dell'onda
%plot(x,u(:,1),'-b',x,u(:,round(Nt/10)),'--g',x,u(:,round(Nt/3)),':b',x,u(:,Nt),'-.r')
%title('square-wave test within the lax Method')
%xlabel('X')
ylabel('Amplitude(X)')