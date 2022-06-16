%Automatización modelo ATP (Basado en paper "Equivalencing the collector
%system of a large wind power plant" E.Muljadi)
clc
clear all

%=========CONSTANTES========
Un=31.5e3;            %tensión nominal
Sb=100e6;             %potencia base
Zb=Un^2/Sb;           %impedancia base
S=[95 120 150 240 400 500 630];                       %mm2
R=[0.411 0.325 0.267 0.164 0.101 0.080 0.06288];      %ohm/km
X=[0.139 0.133 0.128 0.118 0.115 0.111 0.107 0.103];  %ohm/km
C=[0.169 0.185 0.195 0.235 0.25 0.277 0.1 0.345];     %uF/km
f=50;
w=2*pi*f;
B=w.*C*1e-6;    %S/km

%===============================ENTRADAS=================================
%vector desde-hacia
from=[1 2 3 4 5 6 7 8 9 10 11 12]; %numerar cada nodo
to  =[0 1 2 3 3 2 4 6 8 9 9 8];    %como se conectan los nodos (0 es PC)
l   =[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.10 1.20]; %vector de largo, debe coincidir con desde-hacia
sec=[120 120 120 120 120 120 120 120 120 120 120 120]; %vector de secciones, debe coincidir con desde-hacia, secciones permitidas: [95 120 150 240 400 500 630]; 

%================================PROGRAMA=================================
%Se construyen matrices de adyacencia, largo y secciones
from=from+1;
to=to+1;
AD=sparse(from', to', ones((length(to)),1), length(from)+1, length(from)+1);
SEC=sparse(from', to', sec, length(from)+1, length(from)+1);
L=sparse(from', to', l, length(from)+1, length(from)+1);
AD=full(AD);
SEC=full(SEC);
L=full(L);
%visualizacion layuout
names = cell(1,length(AD));
names{1,1}='PC';
for i=2:(length(AD));
    names{1,i}=num2str(i-1);
end

bg = biograph(AD,names,'ShowArrows','off','ShowWeights','off');
h = view(bg);  

%=======CÁLCULO DE RESISTENCIA Y REACTANCIA POR CONDUCTOR======
ru=zeros(length(AD));
xu=zeros(length(AD));
bu=zeros(length(AD));
ultima_bifurcacion=[]; 

for i=1:length(AD)
    for j=1:length(AD)
        if (SEC(i,j)~=0)
        ru(i,j)=R(S==SEC(i,j));
        xu(i,j)=X(S==SEC(i,j)); 
        bu(i,j)=B(S==SEC(i,j));
        end
    end
end

%Construyo matriz con resistencia y reactancia considerando el largo
Ru=ru.*L;
Xu=xu.*L;
Bu=bu.*L;

%===============================================================


%========IDENTIFICAR NODOS CON BIFURCACIONES==============

nodos_bif=[];
sum_AD_nozero=sum(AD ~= 0, 1);   %por cada columna de la matriz, suma los elementos que no son 0
for i=1:length(sum_AD_nozero)    
    if sum_AD_nozero(i)>1        %si la suma > 1, hay mas de una conexión al nodo, hay bifurcación
        nodos_bif=[nodos_bif i-1];
    end
end


%=========================MAIN=================================
%Se va resolviendo aguas abajo hasta no tener ningun nodo con bifurcación
    
while ~isempty(nodos_bif)

sum_AD=sum(AD,1);

%identifico los últimos nodos con bifurcación
ultima_bifurcacion=[];
for i=1:length(nodos_bif) %Para cada nodo con bifurcacion se identifica si es la última
    ultimonodo=0;
    nodos_conectados=[];
    for j=1:length(AD)
        if AD(j,nodos_bif(i)+1)>0 %guardo nodos conectados a la bifurcación
            nodos_conectados=[nodos_conectados j-1];
        end
    end  
            k=1;
            while (k<=length(nodos_conectados)) %Reviso aguas abajo de cada nodo conectado a la bifurcación
                 if (sum_AD(nodos_conectados(k)+1)==0) 
                    ultimonodo=1;               %Si no tienen nada conectado lo identifico
                 else
                     for m=1:length(AD)
                         if AD(m,nodos_conectados(k)+1)>0   %Voy guardando la cadena para veriicar si tiene bifurcación
                            nodos_conectados=[nodos_conectados m-1];
                         end
                     end   
                 end
                 k=k+1;
            end
        
        if (ultimonodo==1) && (isempty(intersect(nodos_conectados, nodos_bif))) %Si algun nodo conectado coincide con los nodos con bifurcación, no es último
            ultima_bifurcacion=[ultima_bifurcacion nodos_bif(i)]; %Si es el último nodo guardo que la bifurcación es la última
        end
end

%Se aplica daisy chain a las últimas bifurcaciones, para ir resolviendo
for i=1: length(ultima_bifurcacion) %Se realiza para todos los nodos que son últimas bifurcaciones
    nodos_conectados=[];
    for j=1:length(AD)
        if AD(j,ultima_bifurcacion(i)+1)>0 %Guardo nodos conectados a la bifurcación
            nodos_conectados=[nodos_conectados j-1];
        end
    end
    
    %Cell con info para conexión en paralelo
    info_circuitos_paralelo=cell(length(nodos_conectados),3);
    
    %Para cada nodo conectado, armo info para daisy chain
    for k=1:length(nodos_conectados)
        aeros_daisy_chain=nodos_conectados(k); %Guardo primer aero conectado
        cantidad_aeros_dc=AD(nodos_conectados(k)+1, ultima_bifurcacion(i)+1);
        R_aeros_daisy_chain=Ru(nodos_conectados(k)+1,ultima_bifurcacion(i)+1);
        X_aeros_daisy_chain=Xu(nodos_conectados(k)+1,ultima_bifurcacion(i)+1);
        ultimonodo=0;
        y=1;
        while (ultimonodo==0) %Hasta no llegar al último aero voy guardando info para resolver daisy chain
            if (sum_AD(aeros_daisy_chain(y)+1)==0)
                ultimonodo=1; %Guardo si es ultimo aero
            else
                for m=1:length(AD)
                    if AD(m,aeros_daisy_chain(y)+1)>0 
                        %Se guarda la info del nodo conectado
                        aeros_daisy_chain=[aeros_daisy_chain m-1];
                        cantidad_aeros_dc=[cantidad_aeros_dc AD(m,aeros_daisy_chain(y)+1)];
                        R_aeros_daisy_chain=[Ru(m,aeros_daisy_chain(y)+1) R_aeros_daisy_chain];
                        X_aeros_daisy_chain=[Xu(m,aeros_daisy_chain(y)+1) X_aeros_daisy_chain];
                        AD(m,aeros_daisy_chain(y)+1)=0; %Se borra info de la matriz de adyacencia (luego se crea columna/fila nueva con info total)
                    end
                end
            end
            y=y+1;
        end
        AD(nodos_conectados(k)+1,ultima_bifurcacion(i)+1)=0; %Se borra info de la matriz de adyacencia
        %============Function daisy chain (segun metodo IID)===============
        [R_dc,X_dc]=daisy_chain(cantidad_aeros_dc, R_aeros_daisy_chain, X_aeros_daisy_chain); %Salida Z equivalente del daisy chain
        %Se guarda la info del daisy para resolver paralelos
        info_circuitos_paralelo{k,1}=sum(cantidad_aeros_dc);
        info_circuitos_paralelo{k,2}=R_dc;
        info_circuitos_paralelo{k,3}=X_dc;
    end
    
    %============Function paralelos (segun metodo IIC)===============
    [R_paralelo, X_paralelo] = circuitos_paralelo(info_circuitos_paralelo); %Impedancia equivalente de todos los circuitos que llegan a la bifurcación
    
    %Se moodifica la matriz de adyacencia y resistencia
    agregar_columna=zeros(length(AD),1); %Se agrega columna con todo ceros ya que nada se conecta aguas abajo
    AD=[AD agregar_columna];
    agregar_fila=zeros(1,length(AD));    %Se agrega fila con un valor N, que corresponde con cuantos aeros se conectan y a que nodo
    total_aeros_bif=0; 
    for q=1:k
        total_aeros_bif=total_aeros_bif+info_circuitos_paralelo{q,1}; %Se cuenta cuantos aeros se conectan al nodo
    end
    agregar_fila(1,ultima_bifurcacion(i)+1)=total_aeros_bif;    %Se agrega a la matriz
    AD=[AD;agregar_fila];
    nodos_bif(nodos_bif==ultima_bifurcacion(i)) = [];           %Se borra de vector el nodo que fue resuelto
    
    %Idem para matriz de resistencias y reactancias
    Ru=[Ru agregar_columna];
    Xu=[Xu agregar_columna];
    agregar_fila(1,ultima_bifurcacion(i)+1)=R_paralelo;         %Se guarda valor devuelto de Function paralelos (segun metodo IIC)
    Ru=[Ru;agregar_fila];
    agregar_fila(1,ultima_bifurcacion(i)+1)=X_paralelo;         %Se guarda valor devuelto de Function paralelos (segun metodo IIC)
    Xu=[Xu;agregar_fila];
end
end

%============TODOS LOS NODOS CON BIFURCACION RESUELTOS===============
%Pendiente resolver el final
 if (ultima_bifurcacion == 1) %Si varios circuitos llegan al PC (indice 1 en matriz), resultado coincide con última salida funcion paralelo
     Rs=R_paralelo/(n^2)/Zb %pu
     Xs=R_paralelo/(n^2)/Zb %pu
 else
     %Queda un último daisy chain por resolver
     %Guardo el nodo conectado al PC
        sum_AD=sum(AD,1);
        nodos_conectados=[];
         for j=1:length(AD)
             if AD(j,1)>0 
                nodos_conectados=[nodos_conectados j-1];
             end
         end
        
        %Guardo info del nodo para resolver DC
        aeros_daisy_chain=nodos_conectados;
        cantidad_aeros_dc=AD(nodos_conectados(1)+1, 1);
        R_aeros_daisy_chain=Ru(nodos_conectados(1)+1,1);
        X_aeros_daisy_chain=Xu(nodos_conectados(1)+1,1);
        ultimonodo=0;
        y=1;
        
        %Voy guardando la cadena
        while (ultimonodo==0) 
            if (sum_AD(aeros_daisy_chain(y)+1)==0) %Si no tiene nodos conectados, último en la cadena
                ultimonodo=1;
            else
                for m=1:length(AD)
                    if AD(m,aeros_daisy_chain(y)+1)>0  %Voy guardando info de la cadena
                        aeros_daisy_chain=[aeros_daisy_chain m-1];
                        cantidad_aeros_dc=[cantidad_aeros_dc AD(m,aeros_daisy_chain(y)+1)];
                        R_aeros_daisy_chain=[Ru(m,aeros_daisy_chain(y)+1) R_aeros_daisy_chain];
                        X_aeros_daisy_chain=[Xu(m,aeros_daisy_chain(y)+1) X_aeros_daisy_chain];
                        
                    end
                end
            end
            y=y+1;
        end
        %============Function daisy chain (segun metodo IID)===============
        %resuelvo el último daisy chain
        [R_dc,X_dc]=daisy_chain(cantidad_aeros_dc, R_aeros_daisy_chain, X_aeros_daisy_chain);
        Rs=R_dc/Zb %pu
        Xs=X_dc/Zb %pu
        
 end
             
     B=sum(Bu(Bu~=0))*Zb %pu             
        

    