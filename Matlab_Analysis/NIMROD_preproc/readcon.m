function [r,z,varr]=readcon(filen,plotit)
%function [r,z,varr]=readcon(filen,plotit)
%reads the contour.bin file produced by nimplot
%input 1 is the number attached to the contour.bin file, i.e. 1000 for
%         contour01000.bin
%input 2 is plotit -> set to 1 to produce contour plots of several 
%      variables for each time slice, otherwise set to zero.
% outputs are 
%    1) r -> 3D array of r locations with indices (ix,iy,ibl)
%    2) z -> 3D array of z locations with indices (ix,iy,ibl)
%    3) varr -> array of varriables with indices 
%        (ix,iy,islice,ibl,variable number)
% where islice is the time slice or poloidal plane number or 
% mode number depending in what is in the contour.bin file
% variable numbers are listed below. To change which variables are
% plotted, alter the line beginning with 'vars='.


file=strcat('contour',num2str(filen,'%05i'),'.bin');
    
fid=fopen(file,'r','b');
bldat=fread(fid,10,'int32');
fclose(fid);
nbl=bldat(2);
nvar=bldat(4);
blx=bldat(7)+1;
bly=bldat(8)+1;
blp=blx*bly;
clear bldat

fid=fopen(file,'r','b');
cdata=fread(fid,inf,'real*4');
fclose(fid);
nslice=(size(cdata,1)-(11+nbl*blp*2+6*(nbl-1)))/(nbl*(blp+2)*nvar);
nslice=max(nslice,1);

for ibl=1:nbl
    r(:,:,ibl)=reshape(cdata(11+(ibl-1)*(2*blp+6):...
        10+(ibl-1)*(2*blp+6)+blp),blx,bly);
    z(:,:,ibl)=reshape(cdata(11+(ibl-1)*(2*blp+6)+blp:...
        10+(ibl-1)*(2*blp+6)+2*blp),blx,bly);
end
cdat2=cdata(11+2*blp*nbl+6*(nbl-1):end);
clear cdata
for sl=1:nslice
    for iv=1:nvar
        for ibl=1:nbl
            %varr(:,:,sl,ibl,iv)=reshape(cdat2(3+(iv-1)*(blp+2)+...
            %    (ibl-1)*(blp+2)*nvar+(sl-1)*nbl*nvar*(blp+2):...
            %    3+(iv-1)*(blp+2)+...
            %(ibl-1)*(blp+2)*nvar+(sl-1)*nbl*nvar*(blp+2)+blp-1),blx,bly);
            varr(:,:,sl,ibl,iv)=reshape(cdat2(3+(iv-1)*(blp+2)+...
                (ibl-1)*(blp+2)*nvar:...
                3+(iv-1)*(blp+2)+...
            (ibl-1)*(blp+2)*nvar+blp-1),blx,bly);
        end
    end
    if (sl<nslice)
    cdat2=cdat2(nbl*nvar*(blp+2)+1:end);
    end
end
clear cdat2
%list of varaible names in the contour.bin file 
labels={...
    'B0_r (T)'          %1
    'B0_z (T)'          %2
    'B0_p_h_i (T)'      %3
    'J0_r (A/m^2)'      %4
    'J0_z (A/m^2)'      %5
    'J0_p_h_i (A/m^2)'  %6
    'V0_r (m/s)'        %7
    'V0_z (m/s)'        %8
    'V0_p_h_i (m/s)'    %9
    'pres0 (Pa)'        %10
    'elec pres0 (Pa)'   %11
    'n0 (m^-^3)'        %12
    'diff shape'        %13
    'B_r (T)'           %14
    'B_z (T)'           %15
    'B_p_h_i (T)'       %16
    'Im B_r (T)'        %17
    'Im B_z (T)'        %18
    'Im B_p_h_i (T)'    %19
    'J_r (A/m^2)'       %20
    'J_z (A/m^2)'       %21
    'J_p_h_i (A/m^2)'   %22
    'Im J_r (A/m^2)'    %23
    'Im J_z (A/m^2)'    %24
    'Im J_p_h_i (A/m^2)' %25
    'V_r (m/s)'         %26
    'V_z (m/s)'         %27
    'V_p_h_i (m/s)'     %28
    'Im V_r (m/s)'      %29
    'Im V_z (m/s)'      %30
    'Im V_p_h_i (m/s)'  %31
    'pressure (Pa)'     %32
    'Im pressure (Pa)'  %33
    'elec pres (Pa)'    %34
    'Im elec pres (Pa)' %35
    'n (m^-^3)'         %36
    'Im n (m^-^3)'      %37
    'conc'              %38
    'Im conc'           %39
    'T_e (eV)'          %40
    'Im T_e (eV)'       %41
    'T_i (eV)'          %42
    'Im T_i (eV)'       %43
    'P_r_a_d (W/m^3)'   %44
    'Im P_r_a_d (W/m^3)' %45
    'P_i_s_o (W/m^3)'   %46
    'Im P_i_s_o (W/m^3)' %47
    'P_i_o_n (W/m^3)'   %48
    'Im P_i_o_n (W/m^3)' %49
    'n_i_m_p (m^-^3)'   %50
    'Im n_i_m_p (m^-^3)' %51
    'n_i_n_o_n (m^-^3)' %52
    'Im n_i_n_o_n (m^-^3)' %53
    'n_n_e_u_t_r_a_l_s (m^-^3)' %54
    'Im n_n_e_u_t_r_a_l_s (m^-^3)' %55
    'n_i_m_p _n_e_u_t_r_a_l'  %56
    'n_i_m_p _+_1' ...           %57
    };
if (plotit)    
% next line determines what set of 12 variables is plotted from
% the list above
    vars=[32 36 42 14 15 16 20 21 22 26 27 28];
%they are plotted in the 3x4 window in this order    
    order=[1 5 9 2 6 10 3 7 11 4 8 12]; 
for sl=1:nslice
    figure(sl+1)
    for iv=1:length(vars)
    if ~isempty(find(varr(:,:,sl,:,vars(iv))))
    h=subplot(3,4,order(iv));
    if ~(minval(varr(:,:,sl,:,vars(iv)))== ...
    maxval(varr(:,:,sl,:,vars(iv))) )
    for ibl=nbl:-1:1
          contour(r(:,:,ibl),z(:,:,ibl),varr(:,:,sl,ibl,vars(iv)),...
          [minval(varr(:,:,sl,:,vars(iv))):...
          (maxval(varr(:,:,sl,:,vars(iv)))-...
           minval(varr(:,:,sl,:,vars(iv))))/40: ...
           maxval(varr(:,:,sl,:,vars(iv)))]);
        hold on     
    end
    end
    axis([minval(r) maxval(r) minval(z) ... 
        maxval(z)]);
    axis image
    set(h,'fontsize',10,'fontname','Book Antiqua');
    xlabel('R (m)','fontsize',10,'fontname','Book Antiqua');
    ylabel('Z (m)','fontsize',10,'fontname','Book Antiqua');
    title(labels(vars(iv)),'fontsize',12,'fontname','Book Antiqua');
    end
    end
    set(sl+1,'position',[130+30*(sl-1) 150-10*(sl-1) 750 800]);
end
end
        

