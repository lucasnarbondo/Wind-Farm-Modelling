function [R_paralelo, X_paralelo] = circuitos_paralelo(info_circuitos_paralelo) %segun sección IIC
size_circuitos_paralelos=size(info_circuitos_paralelo); %Guardo cantidad de circuitos en paralelo
Rsaux=0;
Xsaux=0;
total_aeros=0;
for i=1:size_circuitos_paralelos(1) % (1) porque la cantidad se guarda en las filas
    Rsaux=Rsaux+info_circuitos_paralelo{i,1}^2*info_circuitos_paralelo{i,2}; %Resuelvo según formula IIC en {i,2} se guardan resistencias
    Xsaux=Xsaux+info_circuitos_paralelo{i,1}^2*info_circuitos_paralelo{i,3}; %Resuelvo según formula IIC en {i,3} se guardan reactancias
    total_aeros=total_aeros+info_circuitos_paralelo{i,1};
end
R_paralelo=Rsaux/(total_aeros^2);
X_paralelo=Xsaux/(total_aeros^2);