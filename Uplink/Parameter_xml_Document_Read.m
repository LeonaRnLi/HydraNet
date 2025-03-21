function [SF, BW, Fs, Fc, FFT_factor] = Parameter_xml_Document_Read();
%------------------------------------------------------------------------------------%
% This function aims to read lora parameters in matlab
%
% Ruonan - Jul-06 2023   Established
%------------------------------------------------------------------------------------%

%%=====================Read Tx parameters==================%%

    pathname       = pwd ;
    xmlFileProfile = ['SystemPara','.xml'];
    pathname       = [pathname,'/','SystemPara','/'];
    str            = strcat(pathname,xmlFileProfile);
    xmlDoc         = xmlread(str);
    Tpara          = xmlDoc.getElementsByTagName('Tpara');
    RxparaElement  = Tpara.item(0);
    kk = 1;
    SF = str2num(RxparaElement.item(kk).getTextContent());
    kk = kk+2;
    BW = str2num(RxparaElement.item(kk).getTextContent());
    kk = kk+2;
    Fs = BW;%str2num(RxparaElement.item(kk).getTextContent());
    kk = kk+2;
    Fc = str2num(RxparaElement.item(kk).getTextContent());
    kk = kk+2;
    FFT_factor = str2num(RxparaElement.item(kk).getTextContent());