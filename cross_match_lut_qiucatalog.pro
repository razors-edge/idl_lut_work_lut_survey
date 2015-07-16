;+
; :Author: han
;-
pro cross_match_lut_qiucatalog
close,/all
path='/home/han/Data/lut_data/survey_cat/lut_survey_crossed_20150107/'
lut_file = path + 'lutcat201404_25162_single_match.txt'
catapath = "/home/han/Data/lut_data/survey_cat/"
qiucata_file = catapath + 'nlp_60.cat'
output_path = "/home/han/Data/lut_data/survey_cat/cat_20150107_matchfile/"
match_file = output_path + 'tycho_match.txt'
openw,80,match_file,width = 3000
printf,80,'   lut_ra   lut_dec   lut_sn    lut_m2.0   lut_m2.0err   nlp_ra  nlp_dec   nlp_ab    nlp_b   nlp_v    angula_distance(s)'
nomatch_file = output_path + 'tycho_no_match.txt'
openw,85,nomatch_file,width = 3000
printf,85,'   lut_ra   lut_dec   lut_sn    lut_m2.0   lut_m2.0err   nlp_ra  nlp_dec   nlp_ab    nlp_b   nlp_v    angula_distance(s)'
flare_file = output_path + 'lvs.txt'
openw,90,flare_file,width = 3000
printf,90,'   lut_ra   lut_dec   lut_sn    lut_m2.0   lut_m2.0err   nlp_ra  nlp_dec   nlp_ab    nlp_b   nlp_v    angula_distance(s) '
output_file = output_path + 'match'
print,'lut file ',lut_file
READCOL, lut_file,NAME,RA,DEC,SN2FULL,XCENTER,YCENTER,$
MAG2FULL,MAG2FULLERR,JD,DATETIME,DIRIMG, $
format='A,D,D,f,f,f,f,f,f,A,A'

READCOL, qiucata_file, nlp_ra,nlp_dec,nlp_ab,nlp_b,nlp_v,format='D,D,f,f,f'

entry_device=!d.name
!p.multi=[1,1,1]
set_plot,'ps'
device,file=output_file + '.ps',xsize=8,ysize=6,/inches,xoffset=0.1,yoffset=0.1,/Portrait
device,/color
loadct_plot
!p.position=0



lut_name = NAME
lut_ra = double(RA)
lut_dec = double(DEC)
lut_sn20 = SN2FULL
lut_m20 = MAG2FULL
lut_m20err = MAG2FULLERR
lut_xc = XCENTER
lut_yc = YCENTER
lut_jd = JD
lut_datetime = DATETIME
lut_dirimg = DIRIMG
day = make_array(n_elements(lut_dirimg),/string)
for day_num = 0L, n_elements(lut_dirimg)-1 do begin
day_word = strsplit(lut_dirimg[day_num], '_', /EXTRACT)
day[day_num] = day_word[0]
endfor
day_ind = uniq(day, SORT(day))
;for day_i = 0L , n_elements(day_ind) - 1 do begin
;print, day[day_ind[day_i]]
;endfor

img_ind = uniq(lut_dirimg, SORT(lut_dirimg))
for img_i = 0L , n_elements(img_ind) - 1 do begin
;print, lut_dirimg[img_ind[img_i]]

index_obj = where(lut_dirimg eq lut_dirimg[img_ind[img_i]],count_obj)
if count_obj le 0 then begin
      break
endif else begin
lut_name_obj = lut_name[index_obj]
lut_ra_obj = double(lut_ra[index_obj])
lut_dec_obj = double(lut_dec[index_obj])
lut_sn20_obj = lut_sn20[index_obj]
lut_m20_obj = lut_m20[index_obj]
lut_m20err_obj = lut_m20err[index_obj]
lut_xc_obj = lut_xc[index_obj]
lut_yc_obj = lut_yc[index_obj]
lut_jd_obj = lut_jd[index_obj]
lut_datetime_obj = lut_datetime[index_obj]
lut_dirimg_obj = lut_dirimg[index_obj]

match=match_2d(lut_ra_obj,lut_dec_obj,nlp_ra,nlp_dec,0.017,MATCH_DISTANCE=md_obj)
;print,match
index = where((match ne -1),count)
if count ge 1 then begin
match_index = uniq(match[index])
;print,'count ',count,match_index
angular_dis = make_array(count)
for xi = 0L, count-1 do begin
;print,lut_ra_obj[index[xi]],lut_dec_obj[index[xi]],nlp_ra[match[index[xi]]],nlp_dec[match[index[xi]]]
;angular_distance,lut_ra_obj[index[xi]],lut_dec_obj[index[xi]] $
;,nlp_ra[match[index[xi]]],nlp_dec[match[index[xi]]],angular_dis[xi]
angular_distance,lut_ra_obj[index[xi]],lut_dec_obj[index[xi]] $
,nlp_ra[match[index[xi]]],nlp_dec[match[index[xi]]],ang
angular_dis[xi] = ang
if ang le 0.0014 then begin
print,lut_ra_obj[index[xi]],lut_dec_obj[index[xi]],nlp_ra[match[index[xi]]],nlp_dec[match[index[xi]]],angular_dis[xi],count
printf,80,lut_name_obj[index[xi]],' ',lut_ra_obj[index[xi]],lut_dec_obj[index[xi]],lut_sn20_obj[index[xi]],$
lut_m20_obj[index[xi]],lut_m20err_obj[index[xi]],lut_xc_obj[index[xi]],lut_yc_obj[index[xi]],$
lut_jd_obj[index[xi]],' ',lut_datetime_obj[index[xi]],' ',lut_dirimg_obj[index[xi]],$
nlp_ra[match[index[xi]]],nlp_dec[match[index[xi]]],nlp_ab[match[index[xi]]],$
nlp_b[match[index[xi]]],nlp_v[match[index[xi]]]

;match_obj[5,k] = nlp_ra[match[index[k]]]
;match_obj[6,k] = nlp_dec[match[index[k]]]
;match_obj[7,k] = nlp_ab[match[index[k]]]
;match_obj[8,k] = nlp_b[match[index[k]]]
;match_obj[9,k] = nlp_v[match[index[k]]]
;match_obj[10,k] = match_obj[5,k] * 100000 + match_obj[6,k] * 0.001
;match_obj[11,k] = md[index[k]]


endif
endfor
endif
index1 = where((match eq -1),count1) 
if count1 ge 1 then begin
for yi = 0L, count1-1 do begin
print,lut_ra_obj[index1[yi]],lut_dec_obj[index1[yi]],lut_m20_obj[index1[yi]],lut_dirimg_obj[index1[yi]]
endfor
;print,'count1 ',count1
endif
endelse

endfor

print, 'test memory ',lut_name_obj,lut_dirimg_obj

;READCOL, qiucata_file, nlp_ra,nlp_dec,nlp_ab,nlp_b,nlp_v,format='D,D,f,f,f'
;;match=match_2d(usno_ra,usno_dec,nomad_ra,nomad_dec,.0014,MATCH_DISTANCE=md)

;match=match_2d(lut_ra*3600.,lut_dec*3600.,nlp_ra*3600.,nlp_dec*3600.,20.0,MATCH_DISTANCE=md)
;help,match,md
;index = where((match ne -1),count) 
;print,count
;match_index = uniq(match[index])
;help,match_index

;index1 = where((match eq -1),count1) 
;print,count1


;entry_device=!d.name
;!p.multi=[1,1,1]
;set_plot,'ps'
;device,file=output_file + '.ps',xsize=8,ysize=6,/inches,xoffset=0.1,yoffset=0.1,/Portrait
;device,/color
;loadct_plot
;!p.position=0

;xrange0 = min(lut_m20)-0.5
;xrange1 = max(lut_m20)+0.5
;yrange0 = min(nlp_ab)-0.5
;yrange1 = max(nlp_ab)+0.5
;plot,lut_m20,nlp_ab,$
;psym=8,thick=0.1,$
;color=0,$
;xstyle=1,$
;ystyle=1,$
;xthick=3,ythick=3,$
;xrange=[7,15],$
;yrange=[7,15],$
;position=[0.1,0.1,0.9,1],/noerase,/nodata,$
;xtitle='Observational magnitude (1.5 FWHM)',$
;XCHARSIZE=1.2,$
;ytitle='LUT Fundamental Catalogue magnitude (Qiu, Wu, Han)',$
;YCHARSIZE=1.2
;;axis,xaxis=1,xstyle=1,XCHARSIZE=0.001,xthick=3

;oplot,[7,15],[7,15],thick=5,color=0
;match_obj = make_array(12,count,/DOUBLE)
;for k=0,count-1 do begin 
;if(md[index[k]] le 3.6) then begin
;lut_name = NAME
;lut_ra = double(RA)
;lut_dec = double(DEC)
;lut_sn20 = SN2FULL
;lut_m20 = MAG2FULL
;lut_m20err = MAG2FULLERR
;lut_xc = XCENTER
;lut_yc = YCENTER
;lut_jd = JD
;lut_datetime = DATETIME
;lut_dirimg = DIRIMG

;match_obj[0,k] = lut_ra[index[k]]
;match_obj[1,k] = lut_dec[index[k]]
;match_obj[2,k] = lut_sn20[index[k]]
;match_obj[3,k] = lut_m20[index[k]]
;match_obj[4,k] = lut_m20err[index[k]]
;match_obj[5,k] = nlp_ra[match[index[k]]]
;match_obj[6,k] = nlp_dec[match[index[k]]]
;match_obj[7,k] = nlp_ab[match[index[k]]]
;match_obj[8,k] = nlp_b[match[index[k]]]
;match_obj[9,k] = nlp_v[match[index[k]]]
;match_obj[10,k] = match_obj[5,k] * 100000 + match_obj[6,k] * 0.001
;match_obj[11,k] = md[index[k]]
;;print,name[index[k]],match_obj[*,k]
;endif
;endfor

;multicolumnsort,match_obj,index=10,revert=1
;match_obj_nlp = match_obj[10,*]
;match_ord = uniq(match_obj_nlp)
;match_num = n_elements(match_ord)
;for k=0,match_num-1 do begin  
;  oplot,[match_obj[3,match_ord[k]],match_obj[3,match_ord[k]]],$
;    [match_obj[7,match_ord[k]]+0.02,match_obj[7,match_ord[k]]+0.02],$
;  psym=4,thick=2,$
;  color=2
;  printf,80,match_obj[0,match_ord[k]],match_obj[1,match_ord[k]],$
;    match_obj[2,match_ord[k]],match_obj[3,match_ord[k]],$
;    match_obj[4,match_ord[k]],match_obj[5,match_ord[k]],$
;    match_obj[6,match_ord[k]],match_obj[7,match_ord[k]],$
;    match_obj[8,match_ord[k]],match_obj[9,match_ord[k]],$
;    match_obj[11,match_ord[k]]
;  if(abs(match_obj[2,match_ord[k]] - match_obj[7,match_ord[k]])) ge 1.0 then begin
;  oplot,[match_obj[2,match_ord[k]],match_obj[2,match_ord[k]]],$
;    [match_obj[7,match_ord[k]+0.02],match_obj[7,match_ord[k]]+0.02],$
;  psym=4,thick=3,$
;  color=1
;  printf,90,match_obj[0,match_ord[k]],match_obj[1,match_ord[k]],$
;    match_obj[2,match_ord[k]],match_obj[3,match_ord[k]],$
;    match_obj[4,match_ord[k]],match_obj[5,match_ord[k]],$
;    match_obj[6,match_ord[k]],match_obj[7,match_ord[k]],$
;    match_obj[8,match_ord[k]],match_obj[9,match_ord[k]],$
;    match_obj[11,match_ord[k]]  
;  endif 
;endfor
;print,match_num

;------------------
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
;cgPS2Raster, output_file + '.ps', /JPEG,/PORTRAIT
close,/all
print,'done'
end
