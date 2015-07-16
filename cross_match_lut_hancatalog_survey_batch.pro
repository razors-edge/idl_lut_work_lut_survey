;+
; :Author: han
;-
pro cross_match_lut_hancatalog_survey_batch
close,/all
path_cat='E:\Data\lut_data\survey_cat\cat\'
list=findfile(path_cat+'*.cat')
file_number=n_elements(list)
print,file_number
path='E:\Data\lut_data\survey_cat\'
hancata_file = path + 'typho2_cata_4_lut_candidate_selection_30s.txt'
match_file = path + 'match.txt'
no_match_file = path + 'no_match.txt'
openw,80,match_file,width = 3000
printf,80,'NAME                    RA              DEC             SNR2           '+$
' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
'nlp_ab           nlp_id      angula_distance(s) '
openw,85,no_match_file,width = 3000
printf,85,'NAME                    RA              DEC             SNR2           '+$
' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
'JD        DATETIME               DIRIMG '
flare_file = path + 'lvs.txt'
openw,90,flare_file,width = 3000
printf,90,'NAME                    RA              DEC             SNR2           '+$
' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
'nlp_ab           nlp_id      angula_distance(s) '
output_file = path + 'LUT_survey_mag'

READCOL, hancata_file, nlp_ra,nlp_dec,nlp_ab,format='D,D,f'

entry_device=!d.name
!p.multi=[1,1,1]
set_plot,'ps'
device,file=output_file + '.ps',xsize=8,ysize=6,/inches,xoffset=0.1,yoffset=0.1,/Portrait
device,/color
loadct_plot
!p.position=0

plot,[0,0],[0,0],$
psym=8,thick=0.1,$
color=0,$
xstyle=1,$
ystyle=1,$
xthick=3,ythick=3,$
xrange=[17,6],$
yrange=[17,6],$
position=[0.1,0.1,0.9,1],/noerase,/nodata,$
xtitle='Observational magnitude (2.0 FWHM)',$
XCHARSIZE=1.2,$
ytitle='LUT Fundamental Catalogue magnitude (Qiu, Wu, Han)',$
YCHARSIZE=1.2
;axis,xaxis=1,xstyle=1,XCHARSIZE=0.001,xthick=3

oplot,[7,15.5],[7,15.5],thick=5,color=0

for i=0L, file_number-1 do begin
lut_file=strcompress(list[i],/remove)
print,lut_file
READCOL, lut_file, NAME_MUL, RA_MUL, DEC_MUL, SNR2_MUL, XCENTER_MUL, YCENTER_MUL, MAG2FULL_MUL, $
  MAG2FULLERR_MUL, JD_MUL, DATETIME_MUL, DIRIMG_MUL, $
  format='A,D,D,X,X,X,X,X,X,X,'+$
  'X,X,X,X,X,X,X,X,X,X,'+$
  'X,X,X,X,X,X,f,X,X,X,'+$
  'X,X,f,f,X,X,X,X,X,X,'+$
  'X,X,f,X,X,X,X,X,X,f,'+$
  'X,X,X,X,D,A,A'
;ord1 = uniq(NAME_MUL)
;;help,ord1
;lut_name = NAME_MUL[ord1]
;lut_ra = double(RA_MUL[ord1])
;lut_dec = double(DEC_MUL[ord1])
;lut_snr2 = SNR2_MUL[ord1]
;lut_x = XCENTER_MUL[ord1]
;lut_y = YCENTER_MUL[ord1]
;lut_m20 = MAG2FULL_MUL[ord1]
;lut_m20err = MAG2FULLERR_MUL[ord1]
;lut_jd = JD_MUL[ord1]
;lut_datetime = DATETIME_MUL[ord1]
;lut_dirimg = DIRIMG_MUL[ord1]

lut_name = NAME_MUL
lut_ra = double(RA_MUL)
lut_dec = double(DEC_MUL)
lut_snr2 = SNR2_MUL
lut_x = XCENTER_MUL
lut_y = YCENTER_MUL
lut_m20 = MAG2FULL_MUL
lut_m20err = MAG2FULLERR_MUL
lut_jd = JD_MUL
lut_datetime = DATETIME_MUL
lut_dirimg = DIRIMG_MUL

;match=match_2d(usno_ra,usno_dec,nomad_ra,nomad_dec,.0014,MATCH_DISTANCE=md)
;match=match_2d(lut_ra*3600.,lut_dec*3600.,nlp_ra*3600.,nlp_dec*3600.,20.0,MATCH_DISTANCE=md)
match=match_2d(lut_ra,lut_dec,nlp_ra,nlp_dec,0.02,MATCH_DISTANCE=md)
;help,match,md
index = where((match ne -1),count) 
print,count
match_index = uniq(match[index])
;help,match_index

index1 = where((match eq -1),count1) 
print,count1

str_obj_emp = make_array(4,count,/string)
match_obj_emp = make_array(12,count,/DOUBLE)
for k=0,count-1 do begin 
if(md[index[k]] le (3.0/3600.)) and ( lut_m20[index[k]] le 15.5 ) and ( lut_snr2[index[k]] ge 5 ) then begin
str_obj_emp[0,k] = lut_name[index[k]]
str_obj_emp[1,k] = lut_datetime[index[k]]
str_obj_emp[2,k] = lut_dirimg[index[k]]
str_obj_emp[3,k] = (strtrim(strtrim(string(long64(nlp_ra[match[index[k]]] * 100000000))),1)) + $
  (strtrim(strtrim(string(nlp_dec[match[index[k]]])),1))
match_obj_emp[0,k] = lut_ra[index[k]]
match_obj_emp[1,k] = lut_dec[index[k]]
match_obj_emp[2,k] = lut_snr2[index[k]]
match_obj_emp[3,k] = lut_x[index[k]]
match_obj_emp[4,k] = lut_y[index[k]]
match_obj_emp[5,k] = lut_m20[index[k]]
match_obj_emp[6,k] = lut_m20err[index[k]]
match_obj_emp[7,k] = lut_jd[index[k]]-2456000.
match_obj_emp[8,k] = nlp_ra[match[index[k]]]
match_obj_emp[9,k] = nlp_dec[match[index[k]]]
match_obj_emp[10,k] = nlp_ab[match[index[k]]]
match_obj_emp[11,k] = md[index[k]]
;printf,80,str_obj[0,k],' ',str_obj[1,k],' ',str_obj[2,k],' ',str_obj[3,k],' ',match_obj[*,k]
;  oplot,[match_obj[5,k],match_obj[5,k]],$
;    [match_obj[5,k],match_obj[5,k]],$
;     psym=4,thick=2,$
;     color=2
endif
endfor

badpnts=0
badpnts = Where( str_obj_emp[0,*] eq '', badcount, $
COMPLEMENT=goodpnts, NCOMPLEMENT=goodcount)
IF goodcount GT 0 THEN str_obj = str_obj_emp[*,goodpnts]
IF goodcount GT 0 THEN match_obj = match_obj_emp[*,goodpnts]
;print,observing_time_night_lst
maxrec1=n_elements(str_obj[0,*])

str_no_obj = make_array(3,count1,/string)
no_match_obj = make_array(8,count1,/DOUBLE)
for k=0,count1-1 do begin 
;if(md[index[k]] le (3.0/3600.)) and ( lut_m20[index[k]] le 15.0 ) then begin
str_no_obj[0,k] = lut_name[index1[k]]
str_no_obj[1,k] = lut_datetime[index1[k]]
str_no_obj[2,k] = lut_dirimg[index1[k]]
no_match_obj[0,k] = lut_ra[index1[k]]
no_match_obj[1,k] = lut_dec[index1[k]]
no_match_obj[2,k] = lut_snr2[index1[k]]
no_match_obj[3,k] = lut_x[index1[k]]
no_match_obj[4,k] = lut_y[index1[k]]
no_match_obj[5,k] = lut_m20[index1[k]]
no_match_obj[6,k] = lut_m20err[index1[k]]
no_match_obj[7,k] = lut_jd[index1[k]]-2456000.
printf,85,str_no_obj[*,k] ,no_match_obj[*,k]
;endif
endfor

;multicolumnsort,match_obj,index=10,revert=1
;for k=0,maxrec1-1 do begin 
;printf,80,str_obj[0,k],' ',str_obj[1,k],' ',str_obj[2,k],' ',str_obj[3,k],' ',match_obj[*,k]
;endfor
;match_obj_nlp = str_obj[3,*]
;match_ord = uniq(match_obj_nlp)
;print,maxrec1
;match_num = n_elements(match_ord)
;print,match_num
for k=0,maxrec1-1 do begin 
  printf,80,str_obj[0,k],$
  match_obj[0,k],match_obj[1,k],$
  match_obj[2,k],match_obj[3,k],$
  match_obj[4,k],match_obj[5,k],$
  match_obj[6,k],match_obj[7,k],' ',$
  str_obj[1,k],' ',str_obj[2,k],' ',$
  match_obj[8,k],match_obj[9,k],$
  match_obj[10,k],$
  str_obj[3,k],' ',match_obj[11,k]      
;  if(abs(match_obj[5,k] - match_obj[10,k])) ge 1.0 then begin
  if(match_obj[10,k] - match_obj[5,k]) ge 1.0 then begin
  oplot,[match_obj[5,k],match_obj[5,k]],$
    [match_obj[10,k],match_obj[10,k]],$
     psym=3,thick=15,$
     color=1
  printf,90,str_obj[0,k],$
    match_obj[0,k],match_obj[1,k],$
    match_obj[2,k],match_obj[3,k],$
    match_obj[4,k],match_obj[5,k],$
    match_obj[6,k],match_obj[7,k],' ',$
    str_obj[1,k],' ',str_obj[2,k],' ',$
    match_obj[8,k],match_obj[9,k],$
    match_obj[10,k],$
    str_obj[3,k],' ',match_obj[11,k]
  endif else begin
      oplot,[match_obj[5,k],match_obj[5,k]],$
    [match_obj[10,k],match_obj[10,k]],$
     psym=3,thick=2,$
     color=2
    endelse
endfor
print,maxrec1
endfor

;res1 = make_array(2,count1,/double)
;for m=0,count1-1 do begin 
;;oplot,[obj_ra[index1[m]],obj_ra[index1[m]]],[obj_dec[index1[m]],obj_dec[index1[m]]],$
;;psym=1,thick=1,$
;;color=0
;res1[0,m] = obj_ra[index1[m]]
;res1[1,m] = obj_dec[index1[m]]
;endfor
;
;match_res=match_2d(res1[0,*],res1[1,*],res[0,*],res[1,*],.0056,MATCH_DISTANCE=md)
;help,match_res
;index_res = where((match_res eq -1),count_res) 
;print,count_res
;
;for n=0,count_res-1 do begin 
;oplot,[res[0,index_res[n]],res[0,index_res[n]]],[res[1,index_res[n]],res[1,index_res[n]]],$
;psym=1,thick=1,$
;color=0
;;print,d2dms(res[0,index_res[n]]/15.0),'  ',d2dms(res[1,index_res[n]])
;print,res[0,index_res[n]],res[1,index_res[n]]
;endfor

close,80
device,/close_file
cgPS2Raster, output_file + '.ps', /JPEG,/PORTRAIT
close,/all
print,'done'
end
