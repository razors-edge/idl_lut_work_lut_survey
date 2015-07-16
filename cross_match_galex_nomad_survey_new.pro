;+
; :Author: han
;-
pro cross_match_galex_nomad_survey_new
close,/all
t0=systime(1)&
nomad_path = '/home/han/Data/lut_data/survey_cat/nomad/'
templist=findfile(nomad_path+'temp_*.txt')
templist_number=n_elements(templist)
print,'old temp ',templist_number,templist
if templist_number gt 1 then begin 
    for ms = 0L, templist_number-1 do begin
      temp_file=strcompress(templist[ms],/remove)
      FILE_DELETE, temp_file
    endfor
endif
path_cat='/home/han/Data/lut_data/survey_cat/lut_survey_crossed_20150107/'
path='/home/han/Data/lut_data/survey_cat/cat_20150107_matchfile/'
match_file = path + 'galex_match_nomad.txt'
no_match_file = path + 'galex_no_match_nomad.txt'

openw,80,match_file,width = 3000
;printf,80,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
;'nlp_ab           nlp_id      angula_distance(s) '
openw,85,no_match_file,width = 3000
;printf,85,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG '
flare_file = path + 'lvs_nomad.txt'
openw,90,flare_file,width = 3000
;printf,90,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
;'nlp_ab           nlp_id      angula_distance(s) '
output_file = path + 'Galex_survey_mag_usno'

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
xrange=[250,350],$
yrange=[40,70],$
position=[0.1,0.1,0.9,1],/noerase,/nodata,$
xtitle='RA',$
XCHARSIZE=1.2,$
ytitle='DEC',$
YCHARSIZE=1.2
;axis,xaxis=1,xstyle=1,XCHARSIZE=0.001,xthick=3

;oplot,[7,15.5],[7,15.5],thick=5,color=0
ab_r = [0]
lut_file="/home/han/Data/lut_data/survey_cat/lut_survey_crossed_20150107/galexcat_single_match.txt"
print,'lut file ',lut_file
READCOL, lut_file,RA,DEC,GLON,GLAT,NUVMAG,NUVMAGERR,format='D,D,D,D,f,f'
galex_ra = double(RA)
galex_dec = double(DEC)
galex_glon = double(GLON)
galex_glat = double(GLAT)
galex_nuvmag = NUVMAG
galex_nuvmagerr = NUVMAGERR
count_nomad_obj = 0
count_nomad_obj_1 = 0
for x = 0L, n_elements(galex_ra)-1 do begin
      print,galex_ra[x],galex_dec[x],galex_nuvmag[x]
      exec_py = 'python /home/han/PythonWorkspace/catalogue_local/nomad_mysql_template_extract.py '$
      + nomad_path + ' ' + string(galex_ra[x])+' '+string(galex_dec[x])+'  0.001 -2 19'
      spawn, exec_py
      templist=findfile(nomad_path+'temp_*.txt')
      templist_number=n_elements(templist)
      if templist_number eq 1 then begin 
        temp_file=strcompress(templist[0],/remove)
        print,'temp ',temp_file
        READCOL, temp_file, nomad_ra,nomad_dec,nomad_b,nomad_r,format='D,D,f,X,f'
        if n_elements(nomad_ra) ge 1 then begin
        nomad_br = nomad_b - nomad_r
        nomad_ab = nomad_r + (-1.7371815 + (4.4195631 * nomad_br) - $
                  (1.0513742 * (nomad_br^2.0) ) + (0.098263364 * (nomad_br^ 3.0))) 
        for xd = 0L , n_elements(nomad_ra)-1 do begin
        printf,80,galex_ra[x],$
          galex_dec[x],galex_nuvmag[x],$
          nomad_ra[xd],nomad_dec[xd],$
          nomad_b[xd],nomad_r[xd],$
          nomad_ab[xd],$
          galex_nuvmag[x]-nomad_ab[xd]
        endfor
        count_nomad_obj = count_nomad_obj + 1
	endif else begin       
        count_nomad_obj_1 = count_nomad_obj_1 + 1
        printf,85,galex_ra[x],$
          galex_dec[x],galex_nuvmag[x]
        endelse
        FILE_DELETE, temp_file
      
    endif

percentage = float(x) / float(n_elements(galex_name)) * 100.0
t2=systime(1)
;print,percentage,'% in total is done , ' , $
;(t2-t0) / 60.0 , ' minitus are spent so far'


endfor
;READCOL, qiucata_file, nlp_ra,nlp_dec,nlp_ab,format='D,D,f'
;
;;match=match_2d(lut_ra*3600.,lut_dec*3600.,nlp_ra*3600.,nlp_dec*3600.,20.0,MATCH_DISTANCE=md)
;match=match_2d(lut_ra,lut_dec,nlp_ra,nlp_dec,0.02,MATCH_DISTANCE=md)
;;help,match,md
;index = where((match ne -1),count) 
;print,count
;match_index = uniq(match[index])
;;help,match_index
;
;index1 = where((match eq -1),count1) 
;print,count1
;
;str_obj_emp = make_array(4,count,/string)
;match_obj_emp = make_array(12,count,/DOUBLE)
;for k=0,count-1 do begin 
;if(md[index[k]] le (3.0/3600.)) and ( lut_m20[index[k]] le 15.5 ) and ( lut_snr2[index[k]] ge 5 ) then begin
;str_obj_emp[0,k] = lut_name[index[k]]
;str_obj_emp[1,k] = lut_datetime[index[k]]
;str_obj_emp[2,k] = lut_dirimg[index[k]]
;str_obj_emp[3,k] = (strtrim(strtrim(string(long64(nlp_ra[match[index[k]]] * 100000000))),1)) + $
;  (strtrim(strtrim(string(nlp_dec[match[index[k]]])),1))
;match_obj_emp[0,k] = lut_ra[index[k]]
;match_obj_emp[1,k] = lut_dec[index[k]]
;match_obj_emp[2,k] = lut_snr2[index[k]]
;match_obj_emp[3,k] = lut_x[index[k]]
;match_obj_emp[4,k] = lut_y[index[k]]
;match_obj_emp[5,k] = lut_m20[index[k]]
;match_obj_emp[6,k] = lut_m20err[index[k]]
;match_obj_emp[7,k] = lut_jd[index[k]]-2456000.
;match_obj_emp[8,k] = nlp_ra[match[index[k]]]
;match_obj_emp[9,k] = nlp_dec[match[index[k]]]
;match_obj_emp[10,k] = nlp_ab[match[index[k]]]
;match_obj_emp[11,k] = md[index[k]]
;;printf,80,str_obj[0,k],' ',str_obj[1,k],' ',str_obj[2,k],' ',str_obj[3,k],' ',match_obj[*,k]
;;  oplot,[match_obj[5,k],match_obj[5,k]],$
;;    [match_obj[5,k],match_obj[5,k]],$
;;     psym=4,thick=2,$
;;     color=2
;endif
;endfor
;
;badpnts=0
;badpnts = Where( str_obj_emp[0,*] eq '', badcount, $
;COMPLEMENT=goodpnts, NCOMPLEMENT=goodcount)
;IF goodcount GT 0 THEN str_obj = str_obj_emp[*,goodpnts]
;IF goodcount GT 0 THEN match_obj = match_obj_emp[*,goodpnts]
;;print,observing_time_night_lst
;maxrec1=n_elements(str_obj[0,*])
;
;str_no_obj = make_array(3,count1,/string)
;no_match_obj = make_array(8,count1,/DOUBLE)
;for k=0,count1-1 do begin 
;;if(md[index[k]] le (3.0/3600.)) and ( lut_m20[index[k]] le 15.0 ) then begin
;str_no_obj[0,k] = lut_name[index1[k]]
;str_no_obj[1,k] = lut_datetime[index1[k]]
;str_no_obj[2,k] = lut_dirimg[index1[k]]
;no_match_obj[0,k] = lut_ra[index1[k]]
;no_match_obj[1,k] = lut_dec[index1[k]]
;no_match_obj[2,k] = lut_snr2[index1[k]]
;no_match_obj[3,k] = lut_x[index1[k]]
;no_match_obj[4,k] = lut_y[index1[k]]
;no_match_obj[5,k] = lut_m20[index1[k]]
;no_match_obj[6,k] = lut_m20err[index1[k]]
;no_match_obj[7,k] = lut_jd[index1[k]]-2456000.
;printf,85,str_no_obj[*,k] ,no_match_obj[*,k]
;;endif
;endfor
;
;;multicolumnsort,match_obj,index=10,revert=1
;;for k=0,maxrec1-1 do begin 
;;printf,80,str_obj[0,k],' ',str_obj[1,k],' ',str_obj[2,k],' ',str_obj[3,k],' ',match_obj[*,k]
;;endfor
;;match_obj_nlp = str_obj[3,*]
;;match_ord = uniq(match_obj_nlp)
;;print,maxrec1
;;match_num = n_elements(match_ord)
;;print,match_num
;for k=0,maxrec1-1 do begin 
;  printf,80,str_obj[0,k],$
;  match_obj[0,k],match_obj[1,k],$
;  match_obj[2,k],match_obj[3,k],$
;  match_obj[4,k],match_obj[5,k],$
;  match_obj[6,k],match_obj[7,k],' ',$
;  str_obj[1,k],' ',str_obj[2,k],' ',$
;  match_obj[8,k],match_obj[9,k],$
;  match_obj[10,k],$
;  str_obj[3,k],' ',match_obj[11,k]      
;;  if(abs(match_obj[5,k] - match_obj[10,k])) ge 1.0 then begin
;  if(match_obj[10,k] - match_obj[5,k]) ge 1.0 then begin
;  oplot,[match_obj[5,k],match_obj[5,k]],$
;    [match_obj[10,k],match_obj[10,k]],$
;     psym=3,thick=15,$
;     color=1
;  printf,90,str_obj[0,k],$
;    match_obj[0,k],match_obj[1,k],$
;    match_obj[2,k],match_obj[3,k],$
;    match_obj[4,k],match_obj[5,k],$
;    match_obj[6,k],match_obj[7,k],' ',$
;    str_obj[1,k],' ',str_obj[2,k],' ',$
;    match_obj[8,k],match_obj[9,k],$
;    match_obj[10,k],$
;    str_obj[3,k],' ',match_obj[11,k]
;  endif else begin
;      oplot,[match_obj[5,k],match_obj[5,k]],$
;    [match_obj[10,k],match_obj[10,k]],$
;     psym=3,thick=2,$
;     color=2
;    endelse
;endfor
;print,maxrec1
print,mean(ab_r)
close,80
device,/close_file
;cgPS2Raster, output_file + '.ps', /JPEG,/PORTRAIT
close,/all
&$
t3=systime(1)
print,(t3-t0) / 60.0 , ' minitus'
print,'done'
end
