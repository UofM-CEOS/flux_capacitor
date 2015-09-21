pro check_Tsrfc

; July 2015, - reads raw files
; separates rows and plots specific variables to check
; 
;; ************************ Input VARIABLES  (30 columns)**************************
; 0 - flag
; 1 - year
; 2 - day
; 3 - hhmm
; 4 - power
; 5 - Tpanel
; 6 - Th2o
; 7 - std_power
; 8 - std_panel
; 9 - std_water


;
;set the low flow cutoff (l/min)



idir='/Projects/ArcticNet/2015/Data/UW_pCO2/raw/Tsrfc/'
odir='/Projects/ArcticNet/2015/Data/UW_pCO2/raw/Tsrfc/screened_figures/'

; Set directory of bulkTemplate.sav
;bulktemp=idir+'pCO2_amd_Template.sav'

CD, idir
bulklist=FILE_SEARCH('*.dat',COUNT=bulk_recs,/NOSORT)

;start loop to process each file in the bulk directory
for k=0,bulk_recs-1 do begin 

  ifile=bulklist(k)
  
 openr,1,ifile 
 nlines = long(file_lines(ifile))
 data=fltarr(10,nlines)
 readf, 1, data
 fname=''
 fname='Start_doy_'+ string(data[2,0])
 ; restore, bulktemp
  
;  data=READ_ASCII(ifile, TEMPLATE=pCO2_amd_Template, COUNT=n_bulk)

;  plt_rec=where (data.field01(*) eq 'EQU', rec_count)
;  if(rec_count gt 5) then begin
    Tsrfc=fltarr(1,nlines)
    Tsrfc=data[6,*] 
      
; Create plots showing the values

    plot1 = PLOT(Tsrfc[0,*], TITLE="Tsrfc"+'_'+ fname)  ; 

; set the plot up
    plot1.COLOR="red"

    
; set  axis information
    ax1 = plot1.AXES
    ax1[0].TITLE = 'Num'
    ax1[1].TITLE = 'Tsrfc'
    ax1[2].HIDE = 1 ; hide top X axis
    ax1[3].HIDE = 1 ; hide right Y axis 
; 
   
    plot1.save, odir+fname+'_Tsrfc.png'  

close, /all
endfor

end

