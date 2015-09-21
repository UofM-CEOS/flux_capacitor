pro check_pCO2

; July 2015, - reads raw files
; separates rows and plots specific variables to check
; 
;; ************************ Input VARIABLES  (30 columns)**************************
; 0 - year
; 1 - month
; 2 - day
; 3 - hour
; 4 - minute
; 5 - second
; 6 - equ_temp (C)
; 7 - CO2(um/m)
; 8 - H20 (mm/m)
; 9 - licor_temp (C)
; 10 - licor_press (mbar)
; 11 - equ_press (?)
; 12 - H2O_flow (l/min)
; 13 - licor_flow (?)
; 14 - equ_pump
; 15 - vent_flow
; 16 - cond_temp
; 17 - dry_box_temp
; 18 - press (ctd in bar)
; 19 - temp  (ctd in C)
; 20 - cond  
; 21 - O2(%sat)
; 22 - O2(ppm)
; 23 - pH
; 24 - eH
; 25 - temp2 (C)
; 26 - Latitude
; 27 - Longitude
; 28 - SOG
; 29 - COG


;
;set the low flow cutoff (l/min)



idir='/Projects/ArcticNet/2015/Data/UW_pCO2/raw/pCO2/'
odir='/Projects/ArcticNet/2015/Data/UW_pCO2/raw/pCO2/screened_figures/'

; Set directory of bulkTemplate.sav
bulktemp=idir+'pCO2_amd_Template.sav'

CD, idir
bulklist=FILE_SEARCH('*.txt',COUNT=bulk_recs,/NOSORT)

;start loop to process each file in the bulk directory
for k=0,bulk_recs-1 do begin 

  ifile=bulklist(k)
  fname=strmid(bulklist(k),8,13)
  
  restore, bulktemp
  
  data=READ_ASCII(ifile, TEMPLATE=pCO2_amd_Template, COUNT=n_bulk)

  plt_rec=where (data.field01(*) eq 'EQU', rec_count)
  if(rec_count gt 5) then begin
    co2=fltarr(1,rec_count)
    eTmp=fltarr(1,rec_count)
    h2o=fltarr(1,rec_count)
    irga_P=fltarr(1,rec_count)
    eP=fltarr(1,rec_count)
    irga_fl=fltarr(1,rec_count)
    h2o_fl=fltarr(1,rec_count)
    ctd_t=fltarr(1,rec_count)
    
    co2[0,*]=data.field09(plt_rec)
    eTmp[0,*]=data.field05(plt_rec)
    h2o[0,*]=data.field12(plt_rec)
    irga_P[0,*]=data.field14(plt_rec)
    eP[0,*]=data.field15(plt_rec)
    irga_fl[0,*]=data.field17(plt_rec)
    h2o_fl[0,*]=data.field16(plt_rec)
    ctd_t[0,*]=data.field27(plt_rec)
    
; Create plots showing the values

    plot1 = PLOT(co2[0,*],/BUFFER, TITLE="CO2"+'_'+ fname)  ; 
    plot2 = PLOT(eTmp[0,*],/BUFFER, TITLE="equ_T"+'_'+ fname) ; 
    plot3 = PLOT(h2o[0,*],/BUFFER, TITLE="h2o"+'_'+ fname)  ;
    plot4 = PLOT(irga_P[0,*],/BUFFER, TITLE="irga_p"+'_'+ fname) ; 
    plot5 = PLOT(eP[0,*],/BUFFER, TITLE="eP"+'_'+ fname)  ; 
    plot6 = PLOT(irga_fl[0,*],/BUFFER, TITLE="gas_fl"+'_'+ fname) ; 
    plot7 = PLOT(h2o_fl[0,*],/BUFFER, TITLE="h2o_fl"+'_'+ fname)  ;
    plot8 = PLOT(ctd_t[0,*],/BUFFER, TITLE="ctd_T"+'_'+ fname) ; 

; set the plot up
    plot1.COLOR="red"
    plot2.COLOR="red"
    plot3.COLOR="red"
    plot4.COLOR="red"
    plot5.COLOR="red"
    plot6.COLOR="red"
    plot7.COLOR="red"
    plot8.COLOR="red"
    
; set  axis information
    ax1 = plot1.AXES
    ax1[0].TITLE = 'Num'
    ax1[1].TITLE = 'co2'
    ax1[2].HIDE = 1 ; hide top X axis
    ax1[3].HIDE = 1 ; hide right Y axis 
; 
    ax2 = plot2.AXES
    ax2[0].TITLE = 'Num'
    ax2[1].TITLE = 'equ_T'
    ax2[2].HIDE = 1 ; hide top X axis
    ax2[3].HIDE = 1 ; hide right Y axis 
;    
    ax3 = plot3.AXES
    ax3[0].TITLE = 'Num'
    ax3[1].TITLE = 'h2o'
    ax3[2].HIDE = 1 ; hide top X axis
    ax3[3].HIDE = 1 ; hide right Y axis     
; 
    ax4 = plot4.AXES
    ax4[0].TITLE = 'Num'
    ax4[1].TITLE = 'irga_P'
    ax4[2].HIDE = 1 ; hide top X axis
    ax4[3].HIDE = 1 ; hide right Y axis 
;
   ax5 = plot5.AXES
    ax5[0].TITLE = 'Num'
    ax5[1].TITLE = 'eP'
    ax5[2].HIDE = 1 ; hide top X axis
    ax5[3].HIDE = 1 ; hide right Y axis 
; 
    ax6 = plot6.AXES
    ax6[0].TITLE = 'Num'
    ax6[1].TITLE = 'gas_fl'
    ax6[2].HIDE = 1 ; hide top X axis
    ax6[3].HIDE = 1 ; hide right Y axis 
;    
    ax7 = plot7.AXES
    ax7[0].TITLE = 'Num'
    ax7[1].TITLE = 'h2o_fl'
    ax7[2].HIDE = 1 ; hide top X axis
    ax7[3].HIDE = 1 ; hide right Y axis     
; 
    ax8 = plot8.AXES
    ax8[0].TITLE = 'Num'
    ax8[1].TITLE = 'ctd_T'
    ax8[2].HIDE = 1 ; hide top X axis
    ax8[3].HIDE = 1 ; hide right Y axis 
    
    plot1.save, odir+fname+'_co2.png'  
    plot2.save, odir+fname+'_equ_T.png'    
    plot3.save, odir+fname+'_h2o.png'   
    plot4.save, odir+fname+'_irga_P.png'  
    plot5.save, odir+fname+'_eP.png'  
    plot6.save, odir+fname+'_gas_fl.png'    
    plot7.save, odir+fname+'_h2o_fl.png'   
    plot8.save, odir+fname+'_ctd_T.png'

endif
;check values for std2
  std2_rec=where (data.field01(*) eq 'STD2', std2_count)
  if(std2_count gt 0) then begin
    std2=fltarr(1,std2_count)
    std2_fl=fltarr(1,std2_count)
    std2_P=fltarr(1,std2_count)
    
    std2[0,*]=data.field09(std2_rec)
    std2_fl[0,*]=data.field17(std2_rec)
    std2_P[0,*]=data.field14(std2_rec)
    
    plot9 = PLOT(std2[0,*],/BUFFER, TITLE="std2"+'_'+ fname)  ; 
    plot10 = PLOT(std2_fl[0,*],/BUFFER, TITLE="std2_fl"+'_'+ fname) ; 
    plot11 = PLOT(std2_P[0,*],/BUFFER, TITLE="std2_P"+'_'+ fname)  ;

;     set the plot up
    plot9.COLOR="red"
    plot10.COLOR="red"
    plot11.COLOR="red"
    
; set  axis information
    ax9 = plot9.AXES
    ax9[0].TITLE = 'Num'
    ax9[1].TITLE = 'std2'
    ax9[2].HIDE = 1 ; hide top X axis
    ax9[3].HIDE = 1 ; hide right Y axis 
; 
    ax10 = plot10.AXES
    ax10[0].TITLE = 'Num'
    ax10[1].TITLE = 'std2_fl'
    ax10[2].HIDE = 1 ; hide top X axis
    ax10[3].HIDE = 1 ; hide right Y axis 
;    
    ax11 = plot11.AXES
    ax11[0].TITLE = 'Num'
    ax11[1].TITLE = 'std2_P'
    ax11[2].HIDE = 1 ; hide top X axis
    ax11[3].HIDE = 1 ; hide right Y axis 
      
    plot9.save, odir+fname+'_std2.png'  
    plot10.save, odir+fname+'_std2_fl.png'    
    plot11.save, odir+fname+'_std2_P.png'   
 
endif

;check values for std3
  std3_rec=where (data.field01(*) eq 'STD3', std3_count)
  if(std3_count gt 0) then begin
    std3=fltarr(1,std3_count)
    std3_fl=fltarr(1,std3_count)
    std3_P=fltarr(1,std3_count)
    
    std3[0,*]=data.field09(std3_rec)
    std3_fl[0,*]=data.field17(std3_rec)
    std3_P[0,*]=data.field14(std3_rec)
    
    plot12 = PLOT(std3[0,*],/BUFFER, TITLE="std3"+'_'+ fname)  ; 
    plot13 = PLOT(std3_fl[0,*],/BUFFER, TITLE="std3_fl"+'_'+ fname) ; 
    plot14 = PLOT(std3_P[0,*],/BUFFER, TITLE="std3_P"+'_'+ fname)  ;

;     set the plot up
    plot12.COLOR="red"
    plot13.COLOR="red"
    plot14.COLOR="red"
    
; set  axis information
    ax12 = plot12.AXES
    ax12[0].TITLE = 'Num'
    ax12[1].TITLE = 'std3'
    ax12[2].HIDE = 1 ; hide top X axis
    ax12[3].HIDE = 1 ; hide right Y axis 
; 
    ax13 = plot13.AXES
    ax13[0].TITLE = 'Num'
    ax13[1].TITLE = 'std3_fl'
    ax13[2].HIDE = 1 ; hide top X axis
    ax13[3].HIDE = 1 ; hide right Y axis 
;    
    ax14 = plot14.AXES
    ax14[0].TITLE = 'Num'
    ax14[1].TITLE = 'std3_P'
    ax14[2].HIDE = 1 ; hide top X axis
    ax14[3].HIDE = 1 ; hide right Y axis 
      
    plot12.save, odir+fname+'_std3.png'  
    plot13.save, odir+fname+'_std3_fl.png'    
    plot14.save, odir+fname+'_std3_P.png'   
 
endif

;check values for std4
  std4_rec=where (data.field01(*) eq 'STD4', std4_count)
  if(std4_count gt 0) then begin
    std4=fltarr(1,std4_count)
    std4_fl=fltarr(1,std4_count)
    std4_P=fltarr(1,std4_count)
    
    std4[0,*]=data.field09(std4_rec)
    std4_fl[0,*]=data.field17(std4_rec)
    std4_P[0,*]=data.field14(std4_rec)
    
    plot15 = PLOT(std4[0,*],/BUFFER, TITLE="std4"+'_'+ fname)  ; 
    plot16 = PLOT(std4_fl[0,*],/BUFFER, TITLE="std4_fl"+'_'+ fname) ; 
    plot17 = PLOT(std4_P[0,*],/BUFFER, TITLE="std4_P"+'_'+ fname)  ;

;     set the plot up
    plot15.COLOR="red"
    plot16.COLOR="red"
    plot17.COLOR="red"
    
; set  axis information
    ax15 = plot15.AXES
    ax15[0].TITLE = 'Num'
    ax15[1].TITLE = 'std4'
    ax15[2].HIDE = 1 ; hide top X axis
    ax15[3].HIDE = 1 ; hide right Y axis 
; 
    ax16 = plot16.AXES
    ax16[0].TITLE = 'Num'
    ax16[1].TITLE = 'std4_fl'
    ax16[2].HIDE = 1 ; hide top X axis
    ax16[3].HIDE = 1 ; hide right Y axis 
;    
    ax17 = plot17.AXES
    ax17[0].TITLE = 'Num'
    ax17[1].TITLE = 'std4_P'
    ax17[2].HIDE = 1 ; hide top X axis
    ax17[3].HIDE = 1 ; hide right Y axis 
      
    plot15.save, odir+fname+'_std4.png'  
    plot16.save, odir+fname+'_std4_fl.png'    
    plot17.save, odir+fname+'_std4_P.png'   
 
endif
close, /all
endfor

end

