pro lut_survey_catalogue_cross_self
close,/all
t0=systime(1)&
path_cat='/home/han/Data/lut_data/survey_cat/cat_20150714/'

list=findfile(path_cat+'lut*.cat')
file_number=n_elements(list)
print,file_number
for i=0L, file_number-1 do begin
;for i=0L, 1-1 do begin  
lut_file=strcompress(list[i],/remove)
print,'lut_file: ',lut_file
file_persentage = float((i+1) / file_number) * 100

out_path='/home/han/Data/lut_data/survey_cat/lut_survey_crossed_20150714/'
file_name_word = strsplit(lut_file , '.' , /EXTRACT)
file_name_word1 = strsplit(file_name_word[0] , '/' , /EXTRACT)
file_name = file_name_word1[6]
print,'file_name: ',file_name
single_match_file = out_path + file_name + '_single_match.txt'
multi_match_file = out_path + file_name + '_multi_match.txt'
print,single_match_file,multi_match_file
openw,84,single_match_file,width=300
openw,85,multi_match_file,width=300

READCOL, lut_file, NAME_MUL, RA_MUL, DEC_MUL, SNR2_MUL, XCENTER_MUL, YCENTER_MUL, MAG2FULL_MUL, $
  MAG2FULLERR_MUL, JD_MUL, DATETIME_MUL, DIRIMG_MUL, $
  format='A,D,D,X,X,X,X,X,X,X,'+$
  'X,X,X,X,X,X,X,X,X,X,'+$
  'X,X,X,X,X,X,f,X,X,X,'+$
  'X,X,f,f,X,X,X,X,X,X,'+$
  'X,X,f,X,X,X,X,X,X,f,'+$
  'X,X,X,X,D,A,A'
;record={NAME_MUL:'',RA_MUL:0.0D,DEC_MUL:0.0d,SNR2_MUL:0.0d,XCENTER_MUL:0.0d,$
;YCENTER_MUL:0.0d,MAG2FULL_MUL:0.0d,MAG2FULLERR_MUL:0.0d,$
;JD_MUL:0.0D,DATETIME_MUL:'',DIRIMG_MUL:'',ind:0L,field_ind:0L}
maxrec = n_elements(NAME_MUL)
ind =( FindGen( maxrec )  )
;ind1 = FindGen( maxrec )
ind1 = multisort(RA_MUL,DEC_MUL,/L64)
;help, ind1
;print, ind1

RA_int = fix(RA_MUL)
DEC_int = fix(DEC_MUL)
;print , RA_int , DEC_int
RA_int_uni = uniq(RA_int, SORT(RA_int))
DEC_int_uni = uniq(DEC_int, SORT(DEC_int))
;print , n_elements(RA_int_uni) , n_elements(DEC_int_uni)
;field_grid = make_array(5,(n_elements(RA_int_uni) * n_elements(DEC_int_uni)))
;for rai = 0 , n_elements(RA_int_uni) - 1 do begin
;print, RA_int[RA_int_uni[rai]]
;for deci = 0 , n_elements(DEC_int_uni) - 1 do begin
;print, DEC_int[DEC_int_uni[deci]]
;field_ind = ( rai * n_elements(DEC_int_uni) ) + deci
;field_grid[0,field_ind] = RA_int[RA_int_uni[rai]]
;field_grid[1,field_ind] = DEC_int[DEC_int_uni[deci]]
;field_grid[2,field_ind] = field_ind
;field_grid[3,field_ind] = 0
;field_grid[4,field_ind] = 0
;endfor
;endfor
;print,field_grid

RA_int = fix(RA_MUL)
DEC_int = fix(DEC_MUL)
data = make_array(8,maxrec,/double)
for sor_nun = 0L, maxrec - 1 do begin
data[0,sor_nun] = RA_MUL[ind1[sor_nun]]
data[1,sor_nun] = DEC_MUL[ind1[sor_nun]]
data[2,sor_nun] = ind1[sor_nun]
data[3,sor_nun] = sor_nun 
data[4,sor_nun] = RA_int[ind1[sor_nun]]
data[5,sor_nun] = DEC_int[ind1[sor_nun]]
data[6,sor_nun] = 0
data[7,sor_nun] = data[2,sor_nun]
endfor
print,'data loaded'
;print,data

field_grid = make_array(5,(n_elements(RA_int_uni) * n_elements(DEC_int_uni)))
total_field_grid_num = (n_elements(RA_int_uni) * n_elements(DEC_int_uni))
for rai = 0L , n_elements(RA_int_uni) - 1 do begin
;print, RA_int[RA_int_uni[rai]]
for deci = 0L , n_elements(DEC_int_uni) - 1 do begin
field_ind = ( rai * n_elements(DEC_int_uni) ) + deci
grid_persentage = float(field_ind+1) / float(n_elements(RA_int_uni) * n_elements(DEC_int_uni)) * 100
field_grid[0,field_ind] = RA_int[RA_int_uni[rai]]
field_grid[1,field_ind] = DEC_int[DEC_int_uni[deci]]
field_grid[2,field_ind] = field_ind
field_grid[3,field_ind] = 0
group_ind = where(((data[4,*] eq field_grid[0,field_ind] ) and (data[5,*] eq field_grid[1,field_ind] ) and $
data[6,*] eq 0),group_count)
if group_count ge 1 then begin
field_grid[4,field_ind] = group_count
group_grid_mm = where(((data[4,*] eq field_grid[0,field_ind]-1 ) and (data[5,*] eq field_grid[1,field_ind]-1 )),group_count_mm)
group_grid_om = where(((data[4,*] eq field_grid[0,field_ind] ) and (data[5,*] eq field_grid[1,field_ind]-1 )),group_count_om)
group_grid_pm = where(((data[4,*] eq field_grid[0,field_ind]+1 ) and (data[5,*] eq field_grid[1,field_ind]-1 )),group_count_pm)
group_grid_mo = where(((data[4,*] eq field_grid[0,field_ind]-1 ) and (data[5,*] eq field_grid[1,field_ind] )),group_count_mo)
group_grid_po = where(((data[4,*] eq field_grid[0,field_ind]+1 ) and (data[5,*] eq field_grid[1,field_ind] )),group_count_po)
group_grid_mp = where(((data[4,*] eq field_grid[0,field_ind]-1 ) and (data[5,*] eq field_grid[1,field_ind]+1 )),group_count_mp)
group_grid_op = where(((data[4,*] eq field_grid[0,field_ind] ) and (data[5,*] eq field_grid[1,field_ind]+1 )),group_count_op)
group_grid_pp = where(((data[4,*] eq field_grid[0,field_ind]+1 ) and (data[5,*] eq field_grid[1,field_ind]+1 )),group_count_pp)
group_grid_all = group_ind
group_grid_count_all = group_count
if group_count_mm ge 1 then begin
group_grid_all = [group_grid_all, group_grid_mm]
group_grid_count_all = group_grid_count_all + group_count_mm
endif
if group_count_om ge 1 then begin
group_grid_all = [group_grid_all, group_grid_om]
group_grid_count_all = group_grid_count_all + group_count_om
endif
if group_count_pm ge 1 then begin
group_grid_all = [group_grid_all, group_grid_pm]
group_grid_count_all = group_grid_count_all + group_count_pm
endif
if group_count_mo ge 1 then begin
group_grid_all = [group_grid_all, group_grid_mo]
group_grid_count_all = group_grid_count_all + group_count_mo
endif
if group_count_po ge 1 then begin
group_grid_all = [group_grid_all, group_grid_po]
group_grid_count_all = group_grid_count_all + group_count_po
endif
if group_count_mp ge 1 then begin
group_grid_all = [group_grid_all, group_grid_mp]
group_grid_count_all = group_grid_count_all + group_count_mp
endif
if group_count_op ge 1 then begin
group_grid_all = [group_grid_all, group_grid_op]
group_grid_count_all = group_grid_count_all + group_count_op
endif
if group_count_pp ge 1 then begin
group_grid_all = [group_grid_all, group_grid_pp]
group_grid_count_all = group_grid_count_all + group_count_pp
endif


;print,field_grid[*,field_ind],group_grid_all,group_grid_count_all,$
;group_count,group_count_mm,group_count_om,group_count_pm,group_count_mo,$
;group_count_po,group_count_mp,group_count_op,group_count_pp

for group_grid_sing_ind = 0L , group_grid_count_all-1 do begin
if (data[6,group_grid_all[group_grid_sing_ind]] eq 0) then begin 
angular_distance, data[0,group_grid_all[group_grid_sing_ind]],data[1,group_grid_all[group_grid_sing_ind]],$
data[0,group_grid_all],data[1,group_grid_all],angular_dis
group_grid_match_ind = where(((angular_dis[*] le 0.001)  ),$
group_grid_match_count)
if group_grid_match_count gt 1 then begin
data[6,group_grid_all[group_grid_sing_ind]] = 1
data[6,group_grid_all[group_grid_match_ind]] = 1 
data[7,group_grid_all[group_grid_sing_ind]] = data[2,group_grid_all[group_grid_sing_ind]]
data[7,group_grid_all[group_grid_match_ind]] = data[2,group_grid_all[group_grid_sing_ind]]
data[0,group_grid_all[group_grid_match_ind]] = data[0,group_grid_all[group_grid_sing_ind]]
data[1,group_grid_all[group_grid_match_ind]] = data[1,group_grid_all[group_grid_sing_ind]]
endif
;for group_grid_match_ind = group_grid_sing_ind+1 , group_grid_count_all-1 do begin
;angular_distance, data[0,group_grid_all[group_grid_sing_ind]],data[1,group_grid_all[group_grid_sing_ind]],$
;data[0,group_grid_all[group_grid_match_ind]],data[1,group_grid_all[group_grid_match_ind]],angular_dis
;if (angular_dis le 0.001)   then begin
;data[6,group_grid_all[group_grid_sing_ind]] = 1
;data[6,group_grid_all[group_grid_match_ind]] = 1
;data[7,group_grid_all[group_grid_sing_ind]] = data[2,group_grid_all[group_grid_sing_ind]]
;data[7,group_grid_all[group_grid_match_ind]] = data[2,group_grid_all[group_grid_sing_ind]]
;endif
;endfor
endif
endfor

;for group_ind_num = 0L, group_count-1 do begin
;data[6,group_ind[group_ind_num]] = field_ind
;print,field_ind
;endfor
print, group_grid_count_all
endif

print, 'grid survey is finished by ', grid_persentage, '% (', $
strcompress(string(field_ind + 1 ),/remove), ' in ', strcompress(string(total_field_grid_num),/remove) ,$
RA_int[RA_int_uni[rai]], DEC_int[DEC_int_uni[deci]], ')'
endfor
endfor

;print,data
for ss = 0L , maxrec-1 do begin
radd = double(STRMID(NAME_MUL[data[2,ss]], 4, 2)) * 15.0
ramm = double(STRMID(NAME_MUL[data[2,ss]], 6, 2)) * 15.0
rass = double(STRMID(NAME_MUL[data[2,ss]], 8, 2)) * 15.0
rad = radd + (ramm / 60.) + (rass / 3600.)
decdd = double(STRMID(NAME_MUL[data[2,ss]], 11, 2)) 
decmm = double(STRMID(NAME_MUL[data[2,ss]], 13, 2)) 
decss = double(STRMID(NAME_MUL[data[2,ss]], 15, 2)) 
decd = decdd + (decmm / 60.) + (decss / 3600.)
;print, NAME_MUL[data[7,ss]],data[0,ss],data[1,ss],data[2,ss],data[7,ss]
;NAME_MUL[data[2,ss]],data[0,data[7,ss]],data[1,data[7,ss]],data[7,data[7,ss]]
;print, NAME_MUL[data[7,ss]], ' ', RA_MUL[data[7,ss]], DEC_MUL[data[7,ss]], SNR2_MUL[data[7,ss]], $
;XCENTER_MUL[data[7,ss]], YCENTER_MUL[data[7,ss]], MAG2FULL_MUL[data[7,ss]], $
;  MAG2FULLERR_MUL[data[7,ss]], JD_MUL[data[7,ss]], ' ', DATETIME_MUL[data[7,ss]], ' ',DIRIMG_MUL[data[7,ss]]
printf,85,NAME_MUL[data[7,ss]], ' ', RA_MUL[data[7,ss]], DEC_MUL[data[7,ss]], $
;NAME_MUL[data[2,ss]], ' ', RA_MUL[data[2,ss]], DEC_MUL[data[2,ss]],$
SNR2_MUL[data[2,ss]], $
XCENTER_MUL[data[2,ss]], YCENTER_MUL[data[2,ss]], MAG2FULL_MUL[data[2,ss]], $
  MAG2FULLERR_MUL[data[2,ss]], JD_MUL[data[2,ss]], ' ', DATETIME_MUL[data[2,ss]], ' ',DIRIMG_MUL[data[2,ss]]
endfor

data_single_ind = uniq(data[7,*], SORT(data[7,*]))
for data_sing_i = 0L , n_elements(data_single_ind) - 1 do begin
;for data_sing_i = 0L , maxrec -1 do begin
;printf,84,NAME_MUL[data[7,data_single_ind[data_sing_i]]], ' ', RA_MUL[data[7,data_single_ind[data_sing_i]]],$
; DEC_MUL[data[7,data_single_ind[data_sing_i]]], SNR2_MUL[data[7,data_single_ind[data_sing_i]]], $
;XCENTER_MUL[data[7,data_single_ind[data_sing_i]]], YCENTER_MUL[data[7,data_single_ind[data_sing_i]]], $
;MAG2FULL_MUL[data[7,data_single_ind[data_sing_i]]],  MAG2FULLERR_MUL[data[7,data_single_ind[data_sing_i]]], $
;JD_MUL[data[7,data_single_ind[data_sing_i]]], ' ', DATETIME_MUL[data[7,data_single_ind[data_sing_i]]], $
;' ',DIRIMG_MUL[data[7,data_single_ind[data_sing_i]]]
printf,84,NAME_MUL[data[7,data_single_ind[data_sing_i]]], ' ', RA_MUL[data[7,data_single_ind[data_sing_i]]],$
DEC_MUL[data[7,data_single_ind[data_sing_i]]], SNR2_MUL[data[7,data_single_ind[data_sing_i]]], $
XCENTER_MUL[data[7,data_single_ind[data_sing_i]]], YCENTER_MUL[data[7,data_single_ind[data_sing_i]]], $
MAG2FULL_MUL[data[7,data_single_ind[data_sing_i]]],  MAG2FULLERR_MUL[data[7,data_single_ind[data_sing_i]]], $
JD_MUL[data[7,data_single_ind[data_sing_i]]], ' ', DATETIME_MUL[data[7,data_single_ind[data_sing_i]]], $
' ',DIRIMG_MUL[data[7,data_single_ind[data_sing_i]]]
endfor
print,n_elements(data_single_ind)



print,'file survey is finished by ',file_persentage,'%'
close,84
close,85

endfor
&$
t3=systime(1)
print,(t3-t0) / 60.0 , ' minitus'
close,/all
print,'done'
end
