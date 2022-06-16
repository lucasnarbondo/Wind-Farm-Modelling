function [R_dc, X_dc] = daisy_chain(cantidad_aeros_dc, R_aeros_daisy_chain,X_aeros_daisy_chain) %segun seccion IID
    Rsaux=0; %inicializo variables
    Xsaux=0;
    cantidad_aeros_dc=fliplr(cantidad_aeros_dc); %se modifica el orden del vector
for m=1:length(R_aeros_daisy_chain) %Resuelvo según formula IID (se desprecia las impedancias Zip
    Rsaux=Rsaux+(sum(cantidad_aeros_dc(1:m)))^2*(R_aeros_daisy_chain(m));
    Xsaux=Xsaux+(sum(cantidad_aeros_dc(1:m)))^2*(X_aeros_daisy_chain(m));
end
    R_dc=Rsaux/(sum(cantidad_aeros_dc)^2);
    X_dc=Xsaux/(sum(cantidad_aeros_dc)^2);
end