pro LUT_survey_statistics
close,/all
path = 'E:\WORK\LUT\LUT_data\'
inputfile=filepath('wcsurv23456_radec.dat',$
root_dir = path)

if(n_elements(inputfile) eq 0) then $
message,'Argument FILE is underfined'
maxrec_all=FILE_LINES(inputfile)
array = make_array(4,maxrec_all)

openr,50,inputfile
point_lun,50,0
str = ''
readf,50,str
readf,50,array
close,50

entry_device=!d.name
!p.multi=[1,1,1]
set_plot,'ps'
device,file=outputfile + '.ps',xsize=8,ysize=6,/inches,xoffset=0.1,yoffset=0.1,/Portrait
device,/color
loadct_plot
!p.position=0

;ploting map
;glactc,radeg,decdeg,2000,ragaladeg,decgaladeg,1, /degree
map_set,0,180,/aitoff,/horizon,/noerase
map_grid

device,/close_file
cgPS2Raster, area_time_fig+'.ps', /JPEG,/PORTRAIT

close,/all
device,/close_file
print,'Done'
end