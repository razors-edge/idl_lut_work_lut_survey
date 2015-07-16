;+
; :Author: han
;-
pro lut_survey_selection
close,/all
path='/home/han/Data/lut_data/survey_cat/cat_20150107_matchfile/'
match_file = path + 'match_nomad.txt'
selected_file = path + 'selected_match_nomad.txt'
openw,80,selected_file,width = 3000

;entry_device=!d.name
;!p.multi=[1,1,1]
;set_plot,'ps'
;device,file=selected_file + '.ps',xsize=8,ysize=6,/inches,xoffset=0.1,yoffset=0.1,/Portrait
;device,/color
;loadct_plot
;!p.position=0

;plot,[0,0],[0,0],$
;psym=8,thick=0.1,$
;color=0,$
;xstyle=1,$
;ystyle=1,$
;xthick=3,ythick=3,$
;xrange=[250,350],$
;yrange=[40,70],$
;position=[0.1,0.1,0.9,1],/noerase,/nodata,$
;xtitle='RA',$
;XCHARSIZE=1.2,$
;ytitle='DEC',$
;YCHARSIZE=1.2
;axis,xaxis=1,xstyle=1,XCHARSIZE=0.001,xthick=3

;lut_name_obj[index_nomad_obj[zx]],lut_ra_obj[index_nomad_obj[zx]],$
;          lut_dec_obj[index_nomad_obj[zx]],$
;          lut_sn20_obj[index_nomad_obj[zx]],lut_m20_obj[index_nomad_obj[zx]],$
;          nomad_ra[match_obj[index_nomad_obj[zx]]],nomad_dec[match_obj[index_nomad_obj[zx]]],$
;          nomad_b[match_obj[index_nomad_obj[zx]]],nomad_r[match_obj[index_nomad_obj[zx]]],$
;          nomad_ab[match_obj[index_nomad_obj[zx]]],$
;          lut_m20_obj[index_nomad_obj[zx]]-nomad_ab[match_obj[index_nomad_obj[zx]]],md_obj[index_nomad_obj[zx]],$
;          ' ',lut_dirimg_obj[index_nomad_obj[zx]]


READCOL, match_file, lut_name_obj, lut_ra_obj, lut_dec_obj,lut_sn20_obj,lut_m20_obj,nomad_ra,$
nomad_dec,nomad_b,nomad_r,nomad_ab,mag_dif,md_obj,lut_dirimg_obj, $
  format='A,D,D,f,f,D,D,f,f,f,f,f,A'

index_obj = where(lut_m20_obj lt 10.0,count_obj)
    if count_obj gt 0 then begin
      for i = 0L, count_obj-1 do begin
      printf,80,lut_name_obj[index_obj[i]], lut_ra_obj[index_obj[i]], lut_dec_obj[index_obj[i]],$
      lut_sn20_obj[index_obj[i]],lut_m20_obj[index_obj[i]],nomad_ra[index_obj[i]],$
      nomad_dec[index_obj[i]],nomad_b[index_obj[i]],nomad_r[index_obj[i]],nomad_ab[index_obj[i]],$
      mag_dif[index_obj[i]],md_obj[index_obj[i]],' ',lut_dirimg_obj[index_obj[i]]
      endfor
   endif

close,80
;device,/close_file
;cgPS2Raster, output_file + '.ps', /JPEG,/PORTRAIT
close,/all
print,'done'
end
