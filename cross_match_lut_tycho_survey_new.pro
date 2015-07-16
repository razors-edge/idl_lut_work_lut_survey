;+
; :Author: han
;-
pro cross_match_lut_tycho_survey_new
close,/all
t0=systime(1)&
path='/home/han/Data/lut_data/survey_cat/lut_survey_crossed_20150107/'
lut_file = path + 'lutcat201404_25162_single_match.txt'
path_field='/home/han/Data/lut_data/survey_cat/cat_20150107/wcs_03survey/'
fieldlist=findfile(path_field+'*.txt')
fieldlist_number=n_elements(fieldlist)
print,fieldlist_number
catapath = "/home/han/Data/lut_data/survey_cat/"
qiucata_file = catapath + 'nlp_60.cat'
READCOL, qiucata_file, nlp_ra,nlp_dec,nlp_ab,nlp_b,nlp_v,format='D,D,f,f,f'

path='/home/han/Data/lut_data/survey_cat/cat_20150107_matchfile/'
match_file = path + 'match_tycho.txt'
no_match_file = path + 'no_match_tycho.txt'

openw,80,match_file,width = 3000
;printf,80,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
;'nlp_ab           nlp_id      angula_distance(s) '
openw,85,no_match_file,width = 3000
;printf,85,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG '
flare_file = path + 'lvs_tycho.txt'
openw,90,flare_file,width = 3000
;printf,90,'NAME                    RA              DEC             SNR2           '+$
;' XCENTER         YCENTER         MAG2FULL       MAG2FULLERR   ' + $
;'JD        DATETIME               DIRIMG                  nlp_ra          nlp_dec          '+$
;'nlp_ab           nlp_id      angula_distance(s) '
output_file = path + 'LUT_survey_mag_tycho'

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

READCOL, lut_file,NAME,RA,DEC,SN2FULL,XCENTER,YCENTER,$
MAG2FULL,MAG2FULLERR,JD,DATETIME,DIRIMG, $
format='A,D,D,f,f,f,f,f,f,A,A'

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
print,day,lut_dirimg
day_ind = uniq(day, SORT(day))
for day_i = 0L , n_elements(day_ind) - 1 do begin
fieldfilename = path_field+'wcs_'+day[day_ind[day_i]]+'.txt'
index_field = where(fieldlist eq fieldfilename, count_field)
print,day_i,fieldfilename
if count_field le 0 then begin
  break
endif else begin
  READCOL, fieldfilename, imgname,center_ra,center_dec,format='A,D,D'
  imgname_num  = n_elements(imgname)
  for x = 0, imgname_num -1 do begin
    index_obj = where(lut_dirimg eq imgname[x],count_obj)
    if count_obj le 0 then begin
      break
    endif else begin
      lut_name_obj = lut_name[index_obj]
      lut_ra_obj = lut_ra[index_obj]
      lut_dec_obj = lut_dec[index_obj]
      lut_sn20_obj = lut_sn20[index_obj]
      lut_m20_obj = lut_m20[index_obj]
      lut_m20err_obj = lut_m20err[index_obj]
      lut_xc_obj = lut_xc[index_obj]
      lut_yc_obj = lut_yc[index_obj]
      lut_jd_obj = lut_jd[index_obj]
      lut_datetime_obj = lut_datetime[index_obj]
      lut_dirimg_obj = lut_dirimg[index_obj]
      cata_index = where((nlp_dec le (center_dec[x]+1.5)) and (nlp_dec ge (center_dec[x]-1.5)) ,cata_count)
      tycho_ra = nlp_ra[cata_index]
      tycho_dec = nlp_dec[cata_index]
      tycho_ab = nlp_ab[cata_index]
      tycho_b = nlp_b[cata_index]
      tycho_v = nlp_v[cata_index]      
;        oplot,[min(tycho_ra),min(tycho_ra)],$
;          [min(tycho_dec),min(tycho_dec)],$
;           psym=3,thick=50,$
;           color=0
;        oplot,[min(tycho_ra),min(tycho_ra)],$
;          [max(tycho_dec),max(tycho_dec)],$
;           psym=3,thick=50,$
;           color=0
;        oplot,[max(tycho_ra),max(tycho_ra)],$
;          [min(tycho_dec),min(tycho_dec)],$
;           psym=3,thick=50,$
;           color=0 
;        oplot,[max(tycho_ra),max(tycho_ra)],$
;          [max(tycho_dec),max(tycho_dec)],$
;           psym=3,thick=50,$
;           color=0          
        match_obj=match_2d(lut_ra_obj,lut_dec_obj,tycho_ra,tycho_dec,0.0014,MATCH_DISTANCE=md_obj)
        index_tycho_obj = where(((match_obj ne -1) ), count_tycho_obj)
        print,'match count ',count_tycho_obj
        if count_tycho_obj gt 0 then begin
        for zx = 0L, count_tycho_obj-1 do begin
        print,lut_name_obj[index_tycho_obj[zx]],lut_ra_obj[index_tycho_obj[zx]],$
        lut_dec_obj[index_tycho_obj[zx]]
        printf,80,lut_name_obj[index_tycho_obj[zx]],lut_ra_obj[index_tycho_obj[zx]],$
          lut_dec_obj[index_tycho_obj[zx]],$
          lut_sn20_obj[index_tycho_obj[zx]],lut_m20_obj[index_tycho_obj[zx]],$
          tycho_ra[match_obj[index_tycho_obj[zx]]],tycho_dec[match_obj[index_tycho_obj[zx]]],$
          tycho_b[match_obj[index_tycho_obj[zx]]],tycho_v[match_obj[index_tycho_obj[zx]]],$
          tycho_ab[match_obj[index_tycho_obj[zx]]],$
          lut_m20_obj[index_tycho_obj[zx]]-tycho_ab[match_obj[index_tycho_obj[zx]]],md_obj[index_tycho_obj[zx]],$
          ' ',lut_dirimg_obj[index_tycho_obj[zx]]
         endfor
        oplot,lut_ra_obj[index_tycho_obj],$
          lut_dec_obj[index_tycho_obj],$
           psym=3,thick=2,$
           color=1        
        endif
        index_tycho_obj_1 = where((match_obj eq -1),count_tycho_obj_1) 
        print,'no match count ',count_tycho_obj_1 
        if count_tycho_obj_1 gt 0 then begin
        for zy = 0L, count_tycho_obj_1-1 do begin
        print,zy,lut_name_obj[index_tycho_obj_1[zy]]
        printf,85,lut_name_obj[index_tycho_obj_1[zy]],lut_ra_obj[index_tycho_obj_1[zy]],lut_dec_obj[index_tycho_obj_1[zy]],$
          lut_sn20_obj[index_tycho_obj_1[zy]],lut_m20_obj[index_tycho_obj_1[zy]],' ',lut_dirimg_obj[index_tycho_obj_1[zy]]
        endfor          
        oplot,lut_ra_obj[index_tycho_obj_1],$
          lut_dec_obj[index_tycho_obj_1],$
           psym=3,thick=2,$
           color=3 
        endif
    endelse

;percentage_0 =( (i*(n_elements(day_ind)*n_elements(imgname))) + $
;(day_i * n_elements(imgname)) + x )
;percentage_1 = (file_number*n_elements(day_ind)*n_elements(imgname))
;percentage = float(percentage_0) / float(percentage_1) * 100.0
;t2=systime(1)
;print,percentage,'% in total is done ( ', percentage_0 , ' in ' , percentage_1 , ' ) . ', $
;(t2-t0) / 60.0 , ' minitus are spent so far'

  endfor
endelse
;
;endfor
;
endfor
close,80
device,/close_file
;cgPS2Raster, output_file + '.ps', /JPEG,/PORTRAIT
close,/all
&$
t3=systime(1)
print,(t3-t0) / 60.0 , ' minitus'
print,'done'
end
