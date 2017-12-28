function biweight, data


;  The routine biweight provides an estimator of the location and scale
;  of the data set DATA.  The scale uses the Biweight function in the
;  general formula of "A-estimators." This formula is given on page of
;  416 in UREDA (formula 4). The BIWEIGHT scale estimate is returned as
;  the value XSBIWT. The BIWEIGHT function is given by:
  
;                                  u((1-u*u)**2)     abs(u) <= 1
;                         f(u) =
;                                  0                 abs(u) >  1
;  where u is defined by

;                         u = (XDATA(I) - M) / c*MAD  .

;  M, MAD, and c are the median, the median absolute deviation from
;  the median, and the tuning constant respectively. The tuning
;  constant is a parameter which is chosen depending on the sample
;  size and the specific function being used for the scale estimate.
;  (See page 417 in UREDA).  Here we take c = 9.0.

;--- The biweght location is found using the formula:

;                         T = M + (sums)

;    where M is the sample median and sums are as given on page 421 in UREDA
;    the tuning constant c is set to 6.0 for calculation of the location as 
;    reccommended by Tukey ()

c1 = 6.0 & c2=9.0      ;; tune constants

location = median(data)              ; median of the data
scale = median(abs(data-location))   ; median of absolute deviations

IF scale LE 0.0001 THEN return, [location, scale]

u = (data-location)/(c1*scale)
tmp = u^2
ind_tmp = where(tmp LT 1.0, ntmp)
tmp = (1-tmp)^2
s3 = total((data(ind_tmp)-location)*tmp(ind_tmp))
s4 = total(tmp(ind_tmp))

location = location+s3/s4
FOR i=0, 9 DO BEGIN
    u = (data-location)/(scale*c1)
    tmp = u^2
    ind_tmp = where(tmp LT 1.0, ntmp)
    tmp = (1-tmp)^2
    s3 = total((data(ind_tmp)-location)*tmp(ind_tmp))
    s4 = total(tmp(ind_tmp))
    location = location+s3/s4
ENDFOR

u = (data-location)/(scale*c2)
tmp = u^2
ind_tmp = where(tmp LT 1.0)
s1 = total((data(ind_tmp)-location)^2*(1-tmp(ind_tmp))^4)
s2 = total((1-tmp(ind_tmp))*(1-5*tmp(ind_tmp)))

scale = n_elements(data)/abs(s2)*sqrt(s1/(n_elements(data)-1))  

return, [location, scale]

END

function mad, x, window = window, finite = finite
;+
; NAME:
;   MAD
; PURPOSE:
;   To calculate the Median Absolute Deviation of a set of data in
;   order to calculate the RMS of the noise.
;
; CALLING SEQUENCE:
;   sigma = MAD(X)
;
; INPUTS:
;   X -- A data array
;
; KEYWORD PARAMETERS:
;   None
;
; OUTPUTS:
;   Sigma -- The standard deviation of the data.
;
; MODIFICATION HISTORY:
;
;       Mon Oct 4 2004, Adam Leroy <aleroy@astro>
;               Altered MAD to consider only finite values if
;               the finite keyword is on. NB: you may not always want this.
;
;       Tue Oct 7 15:59:16 2003, Erik Rosolowsky <eros@cosmic>
;		Added Compatibility for distributions with non-zero
;		mean (oops).
;
;       Mon Oct 6 13:26:11 2003, Erik Rosolowsky <eros@cosmic>
;		Written.
;
;
;-
;;; 1./0.674491 is the factor to compute the SDEV for a Gaussian


if (n_elements(window) eq 0) then begin
    if (keyword_set(finite)) then begin
        ind = where(finite(x) eq 1) 
        mad = median(abs(x[ind]-median(x[ind])))/0.6745
    endif else $
      mad = median(abs(x-median(x)))/0.6745
endif else begin
    mad = dblarr(n_elements(x))
    mad = median(abs(x-median(x, window)), window)/0.6745
endelse 

  return, mad
end

FUNCTION mean_ks, val, k, sdev=sdev
;; K-sigma mean using MAD to compute a robust dispersion
;; also with a limited number of points

fin = where(finite(val) EQ 1, nfinite)
IF nfinite LT 3 THEN BEGIN
    IF keyword_set(sdev) THEN return, [!values.f_nan,!values.f_nan] $
    ELSE return, !values.f_nan
ENDIF ELSE BEGIN
    m = median(val)
    s = mad(val)
    FOR i=0,4 DO BEGIN
        ind = where(abs(val-m) LT k*s, nind)
        m = median(val[ind]) & s = mad(val[ind])
    ENDFOR
    IF keyword_set(sdev) THEN $
      return, [mean(val[ind]),stddev(val[ind])] $
    ELSE return, mean(val[ind])
ENDELSE

END

PRO LaunchWidgets
COMMON database, info, frame, calib, ima1, ima2, psf, spectra, points


;; Colors
device, pseudo_color=8, decompose=0
device, CURSOR_STANDARD=60  ;; pointer is a hand
wdelete, !d.window  ;; delete new window ...
; Set to green the last color
loadct,3,/silent
tvlct,red,green,blue,/get
info.ncolors=!d.TABLE_SIZE
red[info.ncolors-1]=0B
green[info.ncolors-1]=255B
blue[info.ncolors-1]=0B
tvlct,red,green,blue

;; Screen size
size_screen = round(get_screen_size()*0.975)

;; Compute the magnification factor

info.magnification = min( [size_screen[0]/$
                           (info.xsize1+info.xsize2+info.xsizePlot)/info.sampling, $
                           size_screen[1]/$
                           (info.ysize12+info.ysizePlot+20.)/info.sampling])
print,'magnification ', info.magnification


;;; Start the interactive part (widgets)
info.base = widget_base(title='IRS low-resolution spectra',/column)
info.base1= widget_base(info.base, /row)
info.base2= widget_base(info.base, /row)
info.base3= widget_base(info.base, /row)
info.draw1 = widget_draw(info.base1 $
                         ,xsize=info.sampling*info.magnification*info.xsize1 $
                         ,ysize=info.sampling*info.magnification*info.ysize12 $
                         ,uvalue="DRAW_WINDOW1",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw2 = widget_draw(info.base2 $
                         ,xsize=info.sampling*info.magnification*info.xsize1 $
                         ,ysize=info.sampling*info.magnification*info.ysizePlot $
                         ,uvalue="DRAW_WINDOW2",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw3 = widget_draw(info.base1 $
                         ,xsize=info.sampling*info.magnification*info.xsize2 $
                         ,ysize=info.sampling*info.magnification*info.ysize12 $
                         ,uvalue="DRAW_WINDOW3",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw4 = widget_draw(info.base2 $
                         ,xsize=info.sampling*info.magnification*info.xsize2 $ 
                         ,ysize=info.sampling*info.magnification*info.ysizePlot $
                         ,uvalue="DRAW_WINDOW4",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw5 = widget_draw(info.base1 $
                         ,xsize=info.sampling*info.magnification*info.xsizePlot $
                         ,ysize=info.sampling*info.magnification*info.ysize12 $
                         ,uvalue="DRAW_WINDOW5",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)
info.draw6 = widget_draw(info.base2 $
                         ,xsize=info.sampling*info.magnification*info.xsizePlot $
                         ,ysize=info.sampling*info.magnification*info.ysizePlot $
                         ,uvalue="DRAW_WINDOW6",retain=2 $
                         ,/button_events, keyboard_events=1,/tracking_events)

;;; Buttons
text=strarr(3)
minscale=-5 & maxscale=5
label_a = widget_label(info.base3,value='MIN',xsize=20)
info.text[0] = widget_text(info.base3,xsize=10,/editable, uvalue="SCALE_MIN",$
                      value=strtrim(info.minscale,2))
label_b = widget_label(info.base3,value='MAX',xsize=20)
info.text[1] = widget_text(info.base3,xsize=10,/editable, uvalue="SCALE_MAX",$
                      value=strtrim(info.maxscale,2))
label_c = widget_label(info.base3,value='RA',xsize=40)
info.text[2] = widget_text(info.base3,xsize=20, value='00:00:00.000')
label_d = widget_label(info.base3,value='DEC',xsize=40)
info.text[3] = widget_text(info.base3,xsize=20, value='+00:00:00.00')

widget_control, /realize, info.base,xoffset=0, yoffset=0
widget_control, info.draw1, get_value=win_id1 & info.win_id1=win_id1
widget_control, info.draw2, get_value=win_id2 & info.win_id2=win_id2
widget_control, info.draw3, get_value=win_id3 & info.win_id3=win_id3
widget_control, info.draw4, get_value=win_id4 & info.win_id4=win_id4
widget_control, info.draw5, get_value=win_id5 & info.win_id5=win_id5
widget_control, info.draw6, get_value=win_id6 & info.win_id6=win_id6
xmanager, 'IrsLow', info.base, /NO_BLOCK

;;; Create the spectra structure

n1=info.xsize2 ;; number of 1st order pixels (2nd window)
n2=info.xsize1 ;; number of 2nd order pixels (1st window)

;; A and B are two components in the case of double spectrum
;; A is normally used (n number of spectra, can be 1 or 2)
spectra = {$
          pos: 0.,$        ;; only one position
          posB: 0., $      ;; 2nd position (in case of close sources)
          n: 0, $
          n1: intarr(n1),$ ;; points used
          n2: intarr(n2),$ ;; points used
          raA: '',$
          decA: '',$
          raB: '',$
          decB: '',$
          w1:  fltarr(n1),$
          f1a: fltarr(n1), $
          f1b: fltarr(n1), $
          e1a: fltarr(n1), $
          e1b: fltarr(n1), $
          w2:  fltarr(n2),$
          f2a: fltarr(n2), $
          f2b: fltarr(n2), $
          e2a: fltarr(n2), $
          e2b: fltarr(n2) $
          }

END

PRO CreateInfo, list, ch
COMMON database

;; NB: read the calibration from the data !!

findpro, 'irslow', dir=path0, /noprint

info = {$
       channel:"long", $
       action:"NONE", $
       caldir:path0[0]+'spice/cal/', $
       aper: 0,$
       out: 0, $
       force: 0, $
       base:0L, $
       base1:0L, $
       base2:0L, $
       base3:0L, $
       text:strarr(4), $
       win_id1:0L, $
       win_id2:0L, $
       win_id3:0L, $
       win_id4:0L, $
       win_id5:0L, $
       win_id6:0L, $
       xsize1:78., $   ;; 2nd order
       xsize2:114., $  ;; 1st order
       ysize12:128., $ ;; 2*info.shift+max(frame.shift)+32
       xsizePlot:50.,$
       ysizePlot:50.,$
       magnification:2., $ 
       sampling:3., $ ;; oversampling factor
       minscale:-2., $
       maxscale:2.,  $ 
       x0:0.0, $        ;; central position
       draw1:0L, $
       draw2:0L, $
       draw3:0L, $
       draw4:0L, $
       draw5:0L, $
       draw6:0L, $
       ncolors:0L,$
       cursor:[0.,0.], $
       tmp:0.,$ ;; temporal variable
       cliplim: fltarr(2,20),$
       unclip: 0.,$
       clipn: 0,$
       order1:0,$
       order2:0,$
       order:0,$
       ra:'',$
       dec:'',$
       nord1: 0,$
       nord2: 0,$
       ord1: intarr(n_elements(list)),$
       ord2: intarr(n_elements(list)),$
       X_: FINDGEN(128) # REPLICATE(1, 128), $ ;; array of x positions
       Y_: REPLICATE(1, 128) # FINDGEN(128), $ ;; array of y positions
       order_mask:intarr(128,128),$
       bmask:intarr(128,128),$
       hot:intarr(128,128),$
       rms:fltarr(128,128),$ ;; background rms
       iw1:fltarr(128),$ ; interpolated wavelength 1st order
       iw2:fltarr(128),$ ; interpolated wavelength 2nd order
       pixsize: 5.12,$
       shift: 192.0, $     ;; shift between same order in two positions
       shift12:0.6 $  ;; further shift between two order centers in two 
                       ;; complementary positions ( empirically measured) ..
       
}

;;; case of short-low
IF ch THEN BEGIN
    info.channel="short"
    info.xsize1 = 80
    info.xsize2 = 126 
    info.pixsize=1.839
    info.shift  =79.0
ENDIF


END



PRO ReadData, list
COMMON database

readmode=strarr(n_elements(list))
for i=0,n_elements(list)-1 do readmode[i]=sxpar(headfits(list[i]),'READMODE')
readmode=strtrim(readmode,2)
ind=where(readmode eq 'RAW',nind)
IF nind GT 0 THEN list = list[ind] $
ELSE BEGIN
    print,'There are no useful files in your list '
    retall
ENDELSE

n = n_elements(list)


;IF n LT 10 THEN BEGIN
;    print,'There are only '+strtrim(string(n),2)+' frames'
;    print,'For less than 10 frames, use irs_low !'
;    IF NOT info.force THEN retall
;ENDIF

print,'There are only '+strtrim(string(n),2)+' frames'


frame=replicate({$
      file:'fileName', $
      aorkey:'aorkey', $
      order:0, $
      ra:0., $
      dec:0., $
      ra_ref:0.,$
      dec_ref:0.,$
      pa:0., $
      image:fltarr(128,128), $   ; original image
      imageb:fltarr(128,128), $  ; image after background subtraction
      mask:bytarr(128,128), $    ; mask sources in single images
      badpixmask:bytarr(128,128), $ ; bad pixel mask
      shift:0.},n)
frame.mask[*,*]=1
frame.badpixmask[*,*]=1

ra = fltarr(n) & dec = fltarr(n) & pa = fltarr(n)
calset = strarr(n)
for i=0,n-1 do begin 
  ima = readfits(list[i],h,/silent) 
  frame[i].file = list[i]
  frame[i].aorkey = sxpar(h,'AORKEY')
  frame[i].ra = sxpar(h,'RA_SLT') 
  frame[i].dec=sxpar(h,'DEC_SLT') 
  frame[i].ra_ref = sxpar(h,'RA_REF') 
  frame[i].dec_ref=sxpar(h,'DEC_REF') 
  frame[i].pa = sxpar(h,'PA_SLT') 
  frame[i].image = ima 
  calset[i] = sxpar(h,'CAL_SET')
  ind = where(finite(frame[i].image) EQ 0, nind)
  IF nind GT 0 THEN frame[i].badpixmask[ind]=0
  fovname=sxpar(h,'fovname') 
  S = STRSPLIT(fovname, '_', /EXTRACT) 
  IF s[2] EQ '2nd' THEN frame[i].order=2 ELSE frame[i].order=1
endfor

calset = strtrim(calset,2)
IF max(calset) NE min(calset) THEN BEGIN
    print,'The files have different calibrations. Consider to separate '+$
          'the AORs'
    retall
ENDIF

cset = max(calset)
info.caldir = info.caldir + strmid(cset,0,strlen(cset)-2)+"/"


;;; Check PA
IF abs(max(frame.pa)-min(frame.pa)) GT 3. THEN BEGIN
    print,'Too much difference in the PA:',abs(max(frame.pa)-min(frame.pa)),$
          '. Consider to separate the AORs'
    yes = 'n'
    read,'Do you want to proceede ? [y/n] ',yes
    IF yes EQ 'n' OR yes EQ 'N' THEN  retall
endif

ind1=where(frame.order EQ 1,nind1)
ind2=where(frame.order EQ 2,nind2)


dalpha = (frame.ra-frame[0].ra)/!RADEG
frame.shift =asin(sin(dalpha)*cos(frame.dec/!RADEG)/sin(frame.pa/!RADEG))*!RADEG*3600./info.pixsize

;;; ADDED !!!! Check new coordinates if correct !!!!
;;; Shifts are  defined as positive wrt the minimum shift
frame.shift = frame.shift - min(frame.shift)


IF nind1 GT 0 AND nind2 GT 0 THEN BEGIN
    i1=ind1[0]
    i2=ind2[0]
    print,'shift 1-2 ',sphdist(frame[i1].ra,frame[i1].dec,$
                               frame[i2].ra,frame[i2].dec,/deg)*3600./info.pixsize 
ENDIF


info.nord1 = nind1
info.nord2 = nind2
IF nind1 GT 0 THEN info.ord1[0:nind1-1]=ind1
IF nind2 GT 0 THEN info.ord2[0:nind2-1]=ind2


IF info.channel EQ "long" THEN BEGIN
    info.order_mask[25:56,10:123]= 1
    info.order_mask[63:94,0:77]= 2
    info.order_mask[64:94,115:127]=3
ENDIF ELSE BEGIN
    info.order_mask[0:37,1:126]=1
    info.order_mask[42:77,0:79]=2
    info.order_mask[42:77,105:127]=3    
ENDELSE
ind1=where(info.order_mask EQ 1)
ind2=where(info.order_mask EQ 2)
ind3=where(info.order_mask EQ 3)

frame[*].badpixmask[where(info.order_mask EQ 0)]=0

;;; Define the mask (long low) ... CHECK IT
IF info.channel EQ "long" THEN BEGIN
    info.bmask(*)=256
    info.bmask[26:55,10:123]=0
    info.bmask[25,37:109]=0
    info.bmask[56,60:123]=0
    info.bmask[57,108:123]=0
    info.bmask[63:94,0:77]=0
    info.bmask[64:96,115:127]=0
ENDIF ELSE BEGIN
    info.bmask(*)=256
    info.bmask[3:30,1:16]=0
    info.bmask[3:31,17:25]=0
    info.bmask[3:32,26:50]=0
    info.bmask[3:34,51:75]=0
    info.bmask[3:34,76:90]=0
    info.bmask[4:34,91:117]=0
    info.bmask[5:34,118:126]=0
    info.bmask[44:73,0:22]=0
    info.bmask[44:74,23:79]=0
    info.bmask[46:77,108:126]=0
ENDELSE

;;; redefine the order mask with bmask ...
ind = where(info.bmask NE 0, nind)
info.order_mask[ind]=0



;writefits,'order.fits',info.order_mask
;writefits,'bmask.fits',info.bmask

;;; PSF 
;;; This part contains 4 PSFs (for SL1, SL2, LL1, and LL2) and the
;;; wavelength scaling (unit to arcsec). The PSF and wav. scaling 
;;; have been empirically derived from observations of the 
;;; calibration star HR2194 
;;; xpsf is given in tenths of arcsec

IF info.channel EQ "long" THEN BEGIN

    coeff_ll1=[3.07127,0.132006]
    xpsf_ll1=(findgen(57)-28.)*0.1
    psf_ll1 = [$
              0.0115818,0.0129328,0.00600795,0.00948248,0.0131225,0.0165108,$
              0.0227347,0.0265714, 0.0310602, 0.0445574,0.0592800,0.0766394,$
              0.0753781,0.0777648, 0.0817771, 0.0836075,0.0605007,0.0750984,$
              0.0974903, 0.145791,  0.201916,  0.266191, 0.380845, 0.490338,$
              0.582351, 0.756640,  0.906597,  0.954670, 0.990638, 0.961817,$
              0.893178, 0.829592,  0.674614,  0.508501, 0.356424, 0.275918,$
              0.207209, 0.166131,  0.136309, 0.0888451, 0.153240, 0.149285,$
              0.148900, 0.127596,  0.133883,  0.108551,0.0901136,0.0815308,$
              0.0580304,0.0377453, 0.0289902, 0.0203496,0.0216201,0.0172422,$
              0.0148173,  0.00000,   0.00000]
    
    coeff_ll2=[4.66684,0.0599461]
    xpsf_ll2=(findgen(51)-25.)*0.1
    psf_ll2 = [$
              0.00000,  0.00000,0.0227415,0.0187635,0.0216823,0.0274293,$
              0.0361808,0.0513072,0.0656436,0.0829470, 0.111041, 0.127177,$
              0.145012, 0.160188, 0.173303, 0.162190, 0.151100, 0.182779,$
              0.237648, 0.326166, 0.445045, 0.582023, 0.765055, 0.900208,$
              0.968570, 0.994295, 0.973966, 0.917324, 0.846249, 0.726052,$
              0.562634, 0.438962, 0.319721, 0.228427, 0.147890, 0.122587,$
              0.102381,0.0983983,0.0982167,0.0984223,0.0967540,0.0765603,$
              0.0592154,0.0482789,0.0334914,0.0259302,0.0185674,0.0164043,$
              0.00307746,  0.00000,  0.00000]    
    psf = {$
          c1:coeff_ll1, $
          x1:xpsf_ll1, $
          y1:psf_ll1,  $
          c2:coeff_ll2, $
          x2:xpsf_ll2, $
          y2:psf_ll2}

ENDIF ELSE BEGIN

    coeff_sl1 =[0.511329,0.211599] 
    xpsf_sl1  =(findgen(51)-25.)*0.1
    psf_sl1   =[$
               0.0050847,0.0105367,0.0210735,0.0248952,0.0314740,0.0377058,$
               0.0505238,0.0627873,0.0715812,0.0792611,0.0819414,0.0778417,$
               0.0693609,0.0652158,0.0672305,0.0847710,0.124802 ,0.180216 ,$
               0.272209 ,0.392790 ,0.538006 ,0.683688 ,0.816583 ,0.923249 ,$
               0.990838 ,1.00022  ,0.966092 ,0.877057 ,0.743465 ,0.607266 ,$
               0.460454 ,0.326309 ,0.225219 ,0.149985 ,0.106702 ,0.0851280,$
               0.0770219,0.0857296,0.0957568,0.102899 ,0.105426 ,0.0952786,$
               0.0804128,0.0625610,0.0447632,0.0307864,0.0203423,0.0129989,$
               0.0115446,0.0087628,0.0072633]
    
    coeff_sl2=[1.80495,0.0503671]
    xpsf_sl2=(findgen(43)-21.)*0.1
    psf_sl2   =[$
               0.0086418,0.0172837,0.0225079,0.0325144,0.0446938,0.0570613,$
               0.0673460,0.0810529,0.0908419,0.0986007, 0.102851, 0.116398,$
               0.151437, 0.201494, 0.270288, 0.401954, 0.526189, 0.668211,$
               0.818281, 0.919169, 0.981232, 0.998027, 0.955874, 0.834191,$
               0.730784, 0.617261, 0.483991, 0.374189, 0.300641, 0.254002,$
               0.229265, 0.207988, 0.190897, 0.176500, 0.155299, 0.120399,$
               0.0916296,0.0687154,0.0560544,0.0422125,0.0316686,0.0289158,$
               0.0144579]
    psf = {$
          c1:coeff_sl1, $
          x1:xpsf_sl1, $
          y1:psf_sl1,  $
          c2:coeff_sl2, $
          x2:xpsf_sl2, $
          y2:psf_sl2}   
ENDELSE


;; Normalize PSF inside apertures of:
;;           8 pix     @ 12 micron for SL
;;           7.172 pix @ 27 micron for LL 

IF info.channel EQ "long" THEN BEGIN
    wavref = 27. & width = 7.172/2.*info.pixsize
ENDIF ELSE BEGIN
    wavref = 12. & width = 8./2.*info.pixsize
ENDELSE

x=psf.x1*(psf.c1[0]+psf.c1[1]*wavref) & y=psf.y1
ind = where(abs(x) LE width)
psf_tot = int_tabulated(x[ind],y[ind])
psf.y1=psf.y1/psf_tot

x=psf.x2*(psf.c2[0]+psf.c2[1]*wavref) & y=psf.y2
ind = where(abs(x) LE width)
psf_tot = int_tabulated(x[ind],y[ind])
psf.y2=psf.y2/psf_tot


END


PRO ReadCalib
COMMON database

;;; distortion, wavelenth calibration
IF info.channel EQ "long" THEN ch="b2" ELSE ch="b0"

readcol,info.caldir+ch+'_wavsamp.tbl',o,xc,yc,w,a,x0,y0,x1,y1,x2,y2,x3,y3,$
          format='i,f,f,f,f,f,f,f,f,f,f,f,f',/silent
wave = readfits(info.caldir+ch+'_wavsamp_wave.fits.gz',/silent)
;;; flux calibration
readcol,info.caldir+ch+'_fluxcon.tbl',lref,flux_conv,cal0,cal1,cal2 $
        ,format='x,f,x,x,x,x,x,f,x,f,x,f,x,f',/silent
readcol,info.caldir+ch+'_slitloss_convert.tbl',wapcorr,fapcorr,format='f,f',/silent

;print, info.caldir+ch+'_wavsamp.tbl'
;;; There is no -1 (IRAF convention). I found it by
;;; recomputing the wavsamp.fits from the wavsamp.tbl
dy = abs(y2-y1)

calib = {$
        x_start:0.5*(x1+x2),$  ;; estremities of the wave rectangle
        y_start:0.5*(y1+y2),$
        x_end:0.5*(x0+x3),$
        y_end:0.5*(y0+y3),$
        ind1:where(o EQ 1,nind1), $
        ind2:where(o EQ 2,nind2), $
        w1:w[where(o EQ 1)], $
        w2:w[where(o EQ 2)], $
        dy1:dy[where(o EQ 1)],$  ;; step in wavelength (in fraction of pixels)
        dy2:dy[where(o EQ 2)],$  ;; step in wavelength (in fraction of pixels)
        wave:wave, $
        flux_conv:flux_conv, $
        cal0:cal0, $
        cal1:cal1, $
        cal2:cal2, $
        lref:lref, $
        wapcorr:wapcorr,$
        fapcorr:fapcorr}

;Fix center
IF info.nord1 GT 0 THEN order=1 ELSE order=2
ind=where(info.order_mask EQ order, nind)
x1_ = min(info.x_[ind]) & x2_ = max(info.x_[ind])
ypix = findgen(128)+0.5
IF info.channel EQ "long" THEN ypix = ypix[0:77] ELSE ypix = ypix[0:79]
IF order EQ 1 THEN begin
    xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
    yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
ENDIF ELSE BEGIN
    xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
    yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
ENDELSE
s=sort(yc)
xc=xc[s] & yc=yc[s]
xpix = interpol(xc,yc,ypix)
min60 = min(abs(ypix-60.),imin60)

info.x0 = mean(frame.shift)+info.shift/info.pixsize+xpix(imin60)-x1_
print,'center is: ',info.x0

END


PRO ComputeCoordinates, pos
COMMON database

;;; Central pixel is defined for the two orders (on ima1, ima2)
;;; find it !!
;if info.channel eq "long" then x0=62.4418 else x0=66.709

;; Compute coordinate of the source taking as reference
;; the median coordinates of the frames of the order

maxn = max([info.nord1, info.nord2],imax)
IF imax EQ 0 THEN ind = info.ord1 ELSE ind = info.ord2
ra0 = mean(frame[ind].ra_ref)
dec0 = mean(frame[ind].dec_ref)
pa  = mean(frame[ind].pa)


pa = pa/!radeg
delta = (pos-info.x0)*info.pixsize/3600.
delta = delta/!radeg
ra0 = ra0/!radeg & dec0=dec0/!radeg

dec = asin(cos(delta)*sin(dec0)+sin(delta)*cos(dec0)*cos(pa))
ra  = ra0+asin(sin(pa)*sin(delta)/cos(dec))
dec = dec * !radeg & ra= ra*!radeg
radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc
IF (ideg lt 0 OR imn LT 0 OR xsc LT 0) THEN sign='-' ELSE sign='+'

a3 = strtrim(string(xsec,format='(f6.3)'),2)
IF strlen(a3) EQ 5 THEN a3='0'+a3
d3 = strtrim(string(abs(xsc),format='(f5.2)'),2)
IF strlen(d3) EQ 4 THEN d3='0'+d3

info.ra = strtrim(string(ihr,format='(i2.2)'),2)+':'+$
          strtrim(string(imin,format='(i2.2)'),2)+':'+a3
info.dec = sign+$
           strtrim(string(abs(ideg),format='(i2.2)'),2)+':'+$
           strtrim(string(abs(imn),format='(i2.2)'),2)+':'+d3

widget_control, info.text[2], set_value=strtrim(info.ra,2)
widget_control, info.text[3], set_value=strtrim(info.dec,2)

END

;;; Save spectra ...
PRO SaveSpectra
COMMON database

;if info.channel eq "long" then x0=62.4418 else x0=66.709
IF spectra.n GT 0 THEN BEGIN
    name = spectra.raA+spectra.decA
    IF abs(spectra.pos-info.x0*info.sampling) LT 3. THEN $
      name=name+'_P' ELSE name=name+'_S'
    if info.channel eq "long" then name = name+'_LL.ASC' $
    else  name = name+'_SL.ASC'
    openw,unit,name,/get_lun
    FOR i=0,n_elements(spectra.w1)-1 DO $
      printf,unit,spectra.w1[i],spectra.f1a[i],spectra.e1a[i],spectra.n1[i],$
             format='("1",3(1x,f14.8),1x,i8)'
    FOR i=0,n_elements(spectra.w2)-1 DO $
      printf,unit,spectra.w2[i],spectra.f2a[i],spectra.e2a[i],spectra.n2[i],$
             format='("2",3(1x,f14.8),1x,i8)'
    close,unit & free_lun,unit
ENDIF
IF spectra.n EQ 2 THEN BEGIN
    name = spectra.raB+spectra.decB
    IF abs(spectra.posB-info.x0*info.sampling) LT 3. THEN $
      name=name+'_P' ELSE name=name+'_S'
    if info.channel eq "long" then name = name+'_LL.ASC' $
    else  name = name+'_SL.ASC'
    openw,unit,name,/get_lun
    FOR i=0,n_elements(spectra.w1)-1 DO $
      printf,unit,spectra.w1[i],spectra.f1b[i],spectra.e1b[i],spectra.n1[i],$
             format='("1",3(1x,f14.8),1x,i8)'
    FOR i=0,n_elements(spectra.w2)-1 DO $
      printf,unit,spectra.w2[i],spectra.f2b[i],spectra.e2b[i],spectra.n2[i],$
             format='("2",3(1x,f14.8),1x,i8)'
    close,unit & free_lun,unit
ENDIF
   

END



PRO FitPoints
COMMON database

IF info.channel EQ "long" THEN wref=27. ELSE wref=12.
IF points.order EQ 1 THEN BEGIN
    wref = psf.c1[0]+psf.c1[1]*wref
    w = psf.c1[0]+psf.c1[1]*points.wpix
    xpsf = psf.x1
    ypsf = psf.y1
ENDIF ELSE BEGIN
    wref = psf.c2[0]+psf.c2[1]*wref
    w = psf.c2[0]+psf.c2[1]*points.wpix
    xpsf = psf.x2
    ypsf = psf.y2
ENDELSE
xpsf_ = xpsf * w    ;; xpsf_ is in arcsec
psf_  = ypsf * wref/w  ;; normalized PSF


;; Normalize wrt wavelength direction
IF points.order EQ 1 THEN $
  dy_  = interpol(calib.dy1, calib.w1, points.wpix) $
ELSE dy_  = interpol(calib.dy2, calib.w2, points.wpix)
psf_ = psf_/dy_


xi = points.x
yi = points.y
wi = points.w
mi = points.m
ind = where(finite(wi) AND finite(yi) AND mi, nind)

IF nind GT 0 THEN BEGIN
    xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind]
    IF spectra.n EQ 1 THEN BEGIN
        ypsf = interpol(psf_,xpsf_,xi)
        den = total(ypsf^2 * wi)
        s = total(ypsf*yi * wi)/den 
        es = sqrt(1./den)
        chi2 = (yi-(s>0)*ypsf)*sqrt(wi) 
        print,'min and max chi2 ',min(chi2),max(chi2)
        ind = where((chi2 gt 5. or chi2 LT -3)  AND wi GT 0., nind)
        IF nind GT 0 THEN BEGIN
            print,'bad points ',nind
            print, chi2[ind]
            points.m[ind]=0
            wi[ind]=0.
            den = total(ypsf^2 * wi)
            s = total(ypsf*yi * wi)/den
            es = sqrt(1./den)
        ENDIF
    ENDIF ELSE IF spectra.n EQ 2 THEN BEGIN
        delta = abs(spectra.pos-spectra.posB)*0.5*info.pixsize/info.sampling
        ypsf1 = interpol(psf_,xpsf_-delta,xi)
        ypsf2 = interpol(psf_,xpsf_+delta,xi)
        ;; fit
        A = total(ypsf1^2*wi)
        B = total(ypsf2^2*wi)
        C = total(ypsf1*ypsf2*wi)
        D = total(ypsf1*yi*wi)
        E = total(ypsf2*yi*wi)
        Den = (A*B-C^2)
        s = (D*B-C*E)/Den
        s2 = (A*E-C*D)/Den
        es = sqrt(B/Den)
        es2 = sqrt(A/Den)
        ;; 5-sigma rejection
        chi2 = (yi-(s>0)*ypsf1-(s2>0)*ypsf2)*sqrt(wi)
        ind = where((chi2 lt -3. OR chi2 GT 5) AND wi GT 0., nind)
        IF nind GT 0 THEN BEGIN
            wi[ind]=0.
            A = total(ypsf1^2*wi)
            B = total(ypsf2^2*wi)
            C = total(ypsf1*ypsf2*wi)
            D = total(ypsf1*yi*wi)
            E = total(ypsf2*yi*wi)
            Den = (A*B-C^2)
            s = (D*B-C*E)/Den
            s2 = (A*E-C*D)/Den
            es = sqrt(B/Den)
            es2 = sqrt(A/Den)
        ENDIF
    ENDIF ELSE nindw=0
ENDIF ELSE BEGIN
    s=0 & es=0
    s2=0 & es2=0
ENDELSE


fc = calib.cal0[points.order-1]+calib.cal1[points.order-1]*$
     (points.wpix-calib.lref[points.order-1])+$
     calib.cal2[points.order-1]*(points.wpix-calib.lref[points.order-1])^2
apcorr = interpol(calib.fapcorr,calib.wapcorr,points.wpix)

;; no apcorr means ...
IF info.aper EQ 0 THEN apcorr=apcorr*0.+1.

IF points.order EQ 1 THEN BEGIN
    spectra.f1a[points.j]=s/fc/calib.flux_conv[points.order-1]/apcorr
    spectra.e1a[points.j]=es/fc/calib.flux_conv[points.order-1]/apcorr
    IF spectra.n EQ 2 THEN BEGIN
        spectra.f1b[points.j]=s2/fc/calib.flux_conv[points.order-1]/apcorr
        spectra.e1b[points.j]=es2/fc/calib.flux_conv[points.order-1]/apcorr
    ENDIF        
ENDIF ELSE BEGIN
    spectra.f2a[points.j]=s/fc/calib.flux_conv[points.order-1]/apcorr
    spectra.e2a[points.j]=es/fc/calib.flux_conv[points.order-1]/apcorr
    IF spectra.n EQ 2 THEN BEGIN
        spectra.f2b[points.j]=s2/fc/calib.flux_conv[points.order-1]/apcorr
        spectra.e2b[points.j]=es2/fc/calib.flux_conv[points.order-1]/apcorr
    ENDIF        
ENDELSE

PlotSpectra
PlotPoints,/ylim,dx=dx

END



PRO PlotPoints, ylim=ylim,dx=dx
;; Plot points (and fit) on window #6
COMMON database

wset, info.win_id6
IF info.channel EQ "long" THEN dw=7.172/27.*points.wpix/2.*info.pixsize $
ELSE  dw=8./12.*points.wpix/2.*info.pixsize

;; Fit
IF info.channel EQ "long" THEN wref=27. ELSE wref=12.
IF points.order EQ 1 THEN BEGIN
    wref = psf.c1[0]+psf.c1[1]*wref
    w = psf.c1[0]+psf.c1[1]*points.wpix
    xpsf = psf.x1
    ypsf = psf.y1
ENDIF ELSE BEGIN
    wref = psf.c2[0]+psf.c2[1]*wref
    w = psf.c2[0]+psf.c2[1]*points.wpix
    xpsf = psf.x2
    ypsf = psf.y2
ENDELSE
xpsf_ = xpsf * w    ;; xpsf_ is in arcsec
psf_  = ypsf * wref/w  ;; normalized PSF
IF points.order EQ 1 THEN $
  dy_  = interpol(calib.dy1, calib.w1, points.wpix) $
ELSE dy_  = interpol(calib.dy2, calib.w2, points.wpix)
psf_ = psf_/dy_

IF points.order EQ 1 THEN s=spectra.f1a[points.j] ELSE s=spectra.f2a[points.j]
IF spectra.n EQ 2 THEN $
  IF points.order EQ 1 THEN s2=spectra.f1b[points.j] ELSE s2=spectra.f2b[points.j]
fc = calib.cal0[points.order-1]+calib.cal1[points.order-1]*$
     (points.wpix-calib.lref[points.order-1])+$
     calib.cal2[points.order-1]*(points.wpix-calib.lref[points.order-1])^2
apcorr = interpol(calib.fapcorr,calib.wapcorr,points.wpix)
;; no apcorr means ...
IF info.aper EQ 0 THEN apcorr = apcorr*0.+1.

s=s*fc*calib.flux_conv[points.order-1]*apcorr
IF spectra.n EQ 2 THEN $
  s2=s2*fc*calib.flux_conv[points.order-1]*apcorr

IF NOT keyword_set(ylim) THEN BEGIN
    points.ylim=[min(points.y,/nan)<0,max(points.y,/nan)]
    IF max(psf_*s) GT points.ylim[1] THEN points.ylim[1]=max(psf_*s)
ENDIF

xlim = max(abs([min(points.x),max(points.x)]))

plot, xma=[2,0], yma=[2,0],points.x, points.y $
      , ps=2, symsize=0.2, color=info.ncolors-2,xr=xlim*[-1,1],/xs,/ys,yr=1.1*points.ylim


ind = where(points.w EQ 0 ,nind)
IF nind GT 0 THEN oplot,points.x[ind],points.y[ind],ps=4, color=info.ncolors-1
ind = where(points.m EQ 0 ,nind)
IF nind GT 0 THEN oplot,points.x[ind],points.y[ind],ps=4, color=info.ncolors-2

;; Oplot fit

IF keyword_set(dx) THEN xpsf_=xpsf_+dx 
IF spectra.n EQ 1 THEN $
  oplot, xpsf_,psf_*s, color=info.ncolors-1 $
ELSE IF spectra.n EQ 2 THEN BEGIN
    delta = abs(spectra.pos-spectra.posB)*0.5*info.pixsize/info.sampling
    oplot, xpsf_-delta,psf_*s, color=info.ncolors-2
    oplot, xpsf_+delta,psf_*s2, color=info.ncolors-1
ENDIF




END



PRO SelectPoints, j, order
;; Select points for a PSF fitting
COMMON database


IF spectra.n EQ 1 THEN $
  pos = spectra.pos/info.sampling $
ELSE IF spectra.n EQ 2 THEN $
  pos = (spectra.pos+spectra.posB)*0.5/info.sampling  

IF order EQ 1 THEN BEGIN
    shift1=info.shift/info.pixsize & shift2 = 0.
    ima=ima1 
    IF info.channel EQ "long" THEN BEGIN
        wpix = info.iw1[10:123] 
    ENDIF  ELSE BEGIN
        wpix = info.iw1[1:126]
    ENDELSE
    spectra.n1[*]=0
ENDIF ELSE IF order EQ 2 THEN BEGIN
    shift1=2*info.shift/info.pixsize+info.shift12
    shift2=info.shift/info.pixsize+info.shift12
    ima=ima2
    IF info.channel EQ "long" THEN BEGIN
        wpix = info.iw2[0:77] 
    ENDIF ELSE BEGIN
        wpix = info.iw2[0:79]
    ENDELSE
    spectra.n2[*]=0
ENDIF

;print,'wpix ',wpix[j]

ind=where(info.order_mask EQ order, nind)
;print,'no frames in order ',order,' is ', nind
x1_ = min(info.x_[ind]) & x2_ = max(info.x_[ind])
y1  = min(info.y_[ind]) & y2  = max(info.y_[ind])
ypix = findgen(128)+0.5
IF order EQ 1 THEN BEGIN
    xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
    yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
    IF info.channel EQ "long" THEN ypix = ypix[10:123]
ENDIF ELSE IF order EQ 2 THEN BEGIN
    xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
    yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
    IF info.channel EQ "long" THEN ypix = ypix[0:77] ELSE ypix = ypix[0:79]
ENDIF

s=sort(yc)
xc=xc[s] & yc=yc[s]
xpix = interpol(xc,yc,ypix)
min60 = min(abs(ypix-60.),imin60)
dpix = xpix - xpix[imin60]

IF info.channel EQ "long" THEN dw=7.172/27.*wpix[j]/2.*info.pixsize $
ELSE  dw=8./12.*wpix[j]/2.*info.pixsize
nx = x2_-x1_+1
xx=findgen(nx)


;; select points (add frame_no and x_no fields)
var = info.rms[x1_:x2_,y1+j]^2
hot = info.hot[x1_:x2_,y1+j]
ind=where(hot EQ 1, nind)
wi = var  ;; inverse variance
IF nind GT 0 THEN var[ind]=!values.f_nan
if info.nord1 eq 0 then begin
    xi = xx+frame[info.ord2[0]].shift+shift2 
    yi = frame[info.ord2[0]].imageb[x1_:x2_,y1+j]
    mi = frame[info.ord2[0]].badpixmask[x1_:x2_,y1+j]
    xn = x1_+xx
    fn = replicate(info.ord2[0],nx)
    for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]   
    for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
    for i=1,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
    for i=1,info.nord2-1 do wi=[wi,var] 
    for i=1,info.nord2-1 DO xn=[xn,x1_+xx] 
    for i=1,info.nord2-1 do fn=[fn,replicate(info.ord2[i],nx)]
endif else begin
    xi = xx+frame[info.ord1[0]].shift+shift1  
    yi = frame[info.ord1[0]].imageb[x1_:x2_,y1+j]
    mi = frame[info.ord1[0]].badpixmask[x1_:x2_,y1+j]
    xn = x1_+xx
    fn = replicate(info.ord1[0],nx)
    for i=1,info.nord1-1 do xi=[xi,xx+frame[info.ord1[i]].shift+shift1]  
    for i=1,info.nord1-1 do yi=[yi,frame[info.ord1[i]].imageb[x1_:x2_,y1+j]]   
    for i=1,info.nord1-1 do mi=[mi,frame[info.ord1[i]].badpixmask[x1_:x2_,y1+j]]   
    for i=1,info.nord1-1 do wi=[wi,var] 
    for i=1,info.nord1-1 DO xn=[xn,x1_+xx] 
    for i=1,info.nord1-1 do fn=[fn,replicate(info.ord1[i],nx)]
    if info.nord2 gt 0 then begin
        for i=0,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]    
        for i=0,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
        for i=0,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
        for i=0,info.nord2-1 do wi=[wi,var] 
        for i=0,info.nord2-1 DO xn=[xn,x1_+xx] 
        for i=0,info.nord2-1 do fn=[fn,replicate(info.ord2[i],nx)]
    endif
ENDELSE

;wi = 1./(wi+abs(yi))  ;;; Add positive signal to the variance
wi = 1./wi  ;; inverse variance
ind = where(finite(wi) EQ 0, nind) ;; put to zero 
IF nind GT 0 THEN wi[ind]=0.
s=sort(xi) & xi=xi[s] & yi=yi[s] & wi=wi[s] & mi=mi[s] & fn=fn[s] & xn=xn[s]
;; xi in arcsec relative to the center and distortion corrected

;;; flux density (per arcsec per wav unit)
yi = yi / info.pixsize ;; /1 pix in wavelength direction
wi = wi * info.pixsize^2

;;; relative distance in arcsec
xi=(xi-pos-dpix[j])*info.pixsize  

IF spectra.n EQ 1 THEN $
  ind = where(abs(xi) LT dw, nindw) $
ELSE IF spectra.n EQ 2 THEN $
  ind = where(abs(xi) LT (dw+abs(spectra.pos-spectra.posB)*$
              info.pixsize*0.5/info.sampling), nindw)

xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind] & fn=fn[ind] & xn=xn[ind]

;; Build the structure points, used for plotting and fitting
points = {$
         order: order, $
         wpix: wpix[j], $
         y1: y1, $
         j: j, $
         x: xi, $
         y: yi, $
         w: wi, $
         m: mi, $
         fn: fn, $
         xn: xn, $
         ylim: [min(yi),max(yi)] $
         }

print,'wavelength: ',points.wpix

DisplayImage
PlotSpectra
PlotPoints

END

PRO PlotSpectra
COMMON database

;print,'plotting spectra ...'

ind1 = where(spectra.w1 LT 33)
IF spectra.n EQ 1 THEN BEGIN
    mins = min([(spectra.f1a-spectra.e1a)[ind1],spectra.f2a-spectra.e2a],/nan)
    maxs = max([(spectra.f1a+spectra.e1a)[ind1],spectra.f2a+spectra.e2a],/nan)
ENDIF ELSE IF spectra.n EQ 2 THEN BEGIN
    mins = min([(spectra.f1a-spectra.e1a)[ind1],spectra.f2a-spectra.e2a,$
               (spectra.f1b-spectra.e1b)[ind1],spectra.f2b-spectra.e2b],/nan)
    maxs = max([(spectra.f1a+spectra.e1a)[ind1],spectra.f2a+spectra.e2a,$
               (spectra.f1b+spectra.e1b)[ind1],spectra.f2b+spectra.e2b],/nan)
ENDIF

;print,'min max s ',  mins,maxs
mins = mins>0

wset,info.win_id4
;print, 'total f1 ', total(abs(spectra.f1a),/nan)
IF total(abs(spectra.f1a),/nan) GT 0 THEN BEGIN
    min1=min(spectra.w1,max=max1)
    dw1 = (max1-min1)/(n_elements(spectra.w1)-1.)*.5
    plot,xma=[0,0],yma=[2,0],spectra.w1,spectra.f1a,/ys,/xs $
         ,xr=[min1-dw1,max1+dw1],yr=[mins,maxs]$
         ,col=info.ncolors-2,psym=4,symsize=0.5
    oplot,spectra.w1,spectra.f1a,col=info.ncolors-2
    oplot,spectra.w1,spectra.f1a-spectra.e1a,col=info.ncolors-2,lines=1
    oplot,spectra.w1,spectra.f1a+spectra.e1a,col=info.ncolors-2,lines=1
    IF spectra.n EQ 2 THEN BEGIN
        oplot,spectra.w1,spectra.f1b,col=info.ncolors-1
        oplot,spectra.w1,spectra.f1b-spectra.e1b,col=info.ncolors-1,lines=1
        oplot,spectra.w1,spectra.f1b+spectra.e1b,col=info.ncolors-1,lines=1
    ENDIF
    IF n_elements(points) GT 0 THEN $
      IF points.order EQ 1 THEN BEGIN
        oplot,spectra.w1[points.j]*[1,1],[mins,maxs], col=info.ncolors-1
        wset,info.win_id3
        mag = info.sampling*info.magnification
        ;print,'j is ',info.xsize2,points.j,info.xsize2-points.j
        plots,(info.xsize2-points.j-.5)*[1,1]*mag,$
              [0,info.ysize12]*mag, col=info.ncolors-1,/dev, lines=0
    ENDIF
ENDIF ELSE erase
wset,info.win_id2
;print,'total f2 ', total(abs(spectra.f2a),/nan)
IF total(abs(spectra.f2a),/nan) GT 0 THEN BEGIN
    min2=min(spectra.w2,max=max2)
    dw2 = (max2-min2)/(n_elements(spectra.w2)-1.)*.5
    plot,xma=[0,0],yma=[2,0],spectra.w2,spectra.f2a,/ys,/xs $
         ,xr=[min2-dw2,max2+dw2],yr=[mins,maxs]$
         ,col=info.ncolors-2,psym=4,symsize=0.5
    oplot,spectra.w2,spectra.f2a,col=info.ncolors-2
    oplot,spectra.w2,spectra.f2a-spectra.e2a,col=info.ncolors-2,lines=1
    oplot,spectra.w2,spectra.f2a+spectra.e2a,col=info.ncolors-2,lines=1
    IF spectra.n EQ 2 THEN BEGIN
        oplot,spectra.w2,spectra.f2b,col=info.ncolors-1
        oplot,spectra.w2,spectra.f2b-spectra.e2b,col=info.ncolors-1,lines=1
        oplot,spectra.w2,spectra.f2b+spectra.e2b,col=info.ncolors-1,lines=1
    ENDIF
    IF n_elements(points) GT 0 THEN $ 
      IF points.order EQ 2 THEN BEGIN
        oplot,spectra.w2[points.j]*[1,1],[mins,maxs], col=info.ncolors-1
        wset,info.win_id1
        mag = info.sampling*info.magnification
        plots,(info.xsize1-points.j-.5)*[1,1]*mag,$
              [0,info.ysize12]*mag, col=info.ncolors-1,/dev, lines=0
    ENDIF
ENDIF ELSE erase

END

PRO AlignOrders, p1, p2
COMMON database

;;; fit the positions
p1=p1*info.sampling
p2=p2*info.sampling
nx = info.ysize12*info.sampling
profile1=fltarr(nx)
profile2=fltarr(nx)
xprofile = findgen(nx)
FOR i=0,nx-1 DO profile1[i]=median(ima1[i,30:239])
FOR i=0,nx-1 DO profile2[i]=median(ima2[i,9:179])

;;; Fit a Gaussian profile
;;; check the highest profile
pos1=(p1-5)>0 & pos2=(p1+5)<nx
x = xprofile[pos1:pos2] & y=profile1[pos1:pos2]
a=[y[5],p1,4]
Result = GAUSSFIT( X, Y, A,nterms=3)
z = (x-a[1])/a[2]
p = a[0]*exp(-z^2/2.)
;;; Oplot the fitted parabola and the limits chosen
wset, info.win_id5
oplot,p,x,col=info.ncolors-1
print,'pos, a[1] ',p1,a[1]
p1 = a[1] ;;; fitted position 

pos1=(p2-5)>0 & pos2=(p2+5)<nx
x = xprofile[pos1:pos2] & y=profile2[pos1:pos2]
a=[y[5],p2,4]
Result = GAUSSFIT( X, Y, A,nterms=3)
z = (x-a[1])/a[2]
p = a[0]*exp(-z^2/2.)
;;; Oplot the fitted parabola and the limits chosen
wset, info.win_id5
oplot,p,x,col=info.ncolors-1
print,'pos, a[1] ',p2,a[1]
p2 = a[1] ;;; fitted position 


print,'Shift is ',(p2-p1)/info.sampling

;info.shift12 = info.shift12 - (p2-p1)
info.shift12 = info.shift12-(p2-p1)/info.sampling
ComputeOrder, 2
DisplayImage


END

PRO ExtractSource, pos, wind
COMMON database


print,'wind,pos ',wind,pos

spectra.n = 1
delvarx, points
;;; fit the position
pos=pos*info.sampling
nx = info.ysize12*info.sampling
profile1=fltarr(nx)
profile2=fltarr(nx)
xprofile = findgen(nx)

IF info.channel EQ "long" THEN $
  FOR i=0,nx-1 DO profile1[i]=median(ima1[i,150:239]) $
ELSE FOR i=0,nx-1 DO profile1[i]=median(ima1[i,30:239])

FOR i=0,nx-1 DO profile2[i]=median(ima2[i,9:179])
;;; Fit a Gaussian profile
;;; check the highest profile
CASE wind OF
    0: BEGIN
        maxi = max([profile1[pos],profile2[pos]],imax)
        IF imax EQ 0 THEN profile=profile1 ELSE profile=profile2
    END
    1: profile=profile1
    2: profile=profile2
ENDCASE
pos1=(pos-5)>0 & pos2=(pos+5)<nx
x = xprofile[pos1:pos2] & y=profile[pos1:pos2]
a=[y[5],pos,4]
Result = GAUSSFIT( X, Y, A,nterms=3)
z = (x-a[1])/a[2]
p = a[0]*exp(-z^2/2.)
;;; Oplot the fitted parabola and the limits chosen
wset, info.win_id5
oplot,p,x,col=info.ncolors-1
print,'pos, a[1] ',pos,a[1]
pos = a[1] ;;; fitted position  (- half pixel)


spectra.pos = pos
ComputeCoordinates, pos/info.sampling
spectra.raA = info.ra
spectra.decA = info.dec


;;; Extraction for the two orders
DisplayImage
pos=(a[1])/info.sampling


FOR order=1,2 DO BEGIN

    IF order EQ 1 THEN BEGIN
        shift1=info.shift/info.pixsize & shift2 = 0.
        ima=ima1 
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw1[10:123] 
;            dwpix = info.idw1[10:123] 
        ENDIF  ELSE BEGIN
            wpix = info.iw1[1:126]
;            dwpix = info.idw1[1:126]
        ENDELSE
        spectra.n1[*]=0
        dy_ = interpol(calib.dy1, calib.w1, wpix)
        tot  = total(abs(ima1[pos*info.sampling,*]))
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        shift1=2*info.shift/info.pixsize+info.shift12
        shift2=info.shift/info.pixsize+info.shift12
        ima=ima2
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw2[0:77] 
;            dwpix = info.idw2[0:77] 
        ENDIF ELSE BEGIN
            wpix = info.iw2[0:79]
;            dwpix = info.idw2[0:79]
        ENDELSE
        spectra.n2[*]=0
        dy_ = interpol(calib.dy2, calib.w2, wpix)
        tot = total(abs(ima2[pos*info.sampling,*]),/nan)
    ENDIF
    nx = info.ysize12*info.sampling
    IF info.channel EQ "long" THEN dw=7.172/27.*wpix/2. $
    ELSE  dw=8./12.*wpix/2.
    spectrum = fltarr(n_elements(dw)) 
    espectrum = fltarr(n_elements(dw)) 

    IF tot GT 0 THEN BEGIN
    print,'order tot ',order, tot
        
;;; Plot aperture on the source centered on the profile

        
;;; triple for displaying 
        ypix = findgen(n_elements(dw)*info.sampling)
        dw_= fltarr(n_elements(ypix))
        ii=indgen(n_elements(dw))*info.sampling 
        dw_[ii]=dw & dw_[ii+1]=dw & dw_[ii+2]=dw
        
        IF order EQ 1 THEN wset, info.win_id3 ELSE wset, info.win_id1
    
;;; Compute the flux inside the dw limits by fitting the PSF
;;; to the points ...
        
        indord=where(info.order_mask EQ order, nindord)
        ;print,'no frames in order ',order,' is ', nindord
        x1_ = min(info.x_[indord]) & x2_ = max(info.x_[indord])
        y1  = min(info.y_[indord]) & y2  = max(info.y_[indord])
        ;print,'y1 y2 ',y1,y2
        ypix = findgen(128)+0.5
        IF order EQ 1 THEN BEGIN
            xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
            yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
            IF info.channel EQ "long" THEN ypix = ypix[10:123] ELSE ypix=ypix[1:126]
        ENDIF ELSE IF order EQ 2 THEN BEGIN
            xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
            yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
            IF info.channel EQ "long" THEN ypix = ypix[0:77] ELSE  ypix = ypix[0:79]
        ENDIF
        
        
        s=sort(yc)
        xc=xc[s] & yc=yc[s]
        xpix = interpol(xc,yc,ypix)
        min60 = min(abs(ypix-60.),imin60)
        dpix = xpix - xpix[imin60]
        

        dw = dw * info.pixsize
    
        xx=findgen(x2_-x1_+1)
        x1=fltarr(x2_-x1_+1)+1.
        n1 = n_elements(x1) 

;        ddx = fltarr(n_elements(wpix))*info.magnification*info.sampling
;        wtot=[0]
;        xtot=[0]
;        ytot=[0]
;        stot=[0]
;        ftot=[" "]
        IF info.channel EQ "long" THEN wref=27. ELSE wref=12.
        IF order EQ 1 THEN BEGIN
            wref = psf.c1[0]+psf.c1[1]*wref
            xpsf = psf.x1
            ypsf = psf.y1
        ENDIF ELSE BEGIN
            wref = psf.c2[0]+psf.c2[1]*wref 
            xpsf = psf.x2
            ypsf = psf.y2
        ENDELSE
       
        FOR j=0,y2-y1 DO BEGIN   
            IF j GT 0 THEN BEGIN 
                mag = info.sampling*info.magnification
                plots, col=info.ncolors-1,/dev, lines=0, $
                       ((y2-y1)-[j-1,j])* mag, (pos-dw[j-1:j]/info.pixsize)*mag
                plots, col=info.ncolors-1,/dev, lines=0, $
                       ((y2-y1)-[j-1,j])* mag, (pos+dw[j-1:j]/info.pixsize)*mag
            ENDIF
            IF order EQ 1 THEN BEGIN
                w = psf.c1[0]+psf.c1[1]*wpix[j]
            ENDIF ELSE BEGIN
                w = psf.c2[0]+psf.c2[1]*wpix[j]
            ENDELSE
            xpsf_ = xpsf * w    ;; xpsf_ is in arcsec
            psf_  = ypsf * wref/w  ;; normalized PSF
            psf_ = psf_/dy_[j]  ;; normalize along wavelength direction

            ;; select points
            var = info.rms[x1_:x2_,y1+j]^2
            hot = info.hot[x1_:x2_,y1+j]
            ind=where(hot EQ 1, nind)
            IF nind GT 0 THEN var[ind]=!values.f_nan
            wi = var  ;; inverse variance
            if info.nord1 eq 0 then begin
                xi = xx+frame[info.ord2[0]].shift+shift2 
                yi = frame[info.ord2[0]].imageb[x1_:x2_,y1+j]
                mi = frame[info.ord2[0]].badpixmask[x1_:x2_,y1+j]
                ;si = frame[info.ord2[0]].shift*x1
                ;fi = frame[info.ord2[0]].file
                ;FOR k=1,n1-1 do fi=[fi,frame[info.ord2[0]].file]
                for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]   
                for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
                for i=1,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
                for i=1,info.nord2-1 do wi=[wi,var] 
                ;for i=1,info.nord2-1 do si=[si,frame[info.ord2[i]].shift*x1]
                ;for i=1,info.nord2-1 DO FOR k=0,n1-1 do fi=[fi,frame[info.ord2[i]].file]
            endif else begin
                xi = xx+frame[info.ord1[0]].shift+shift1  
                yi = frame[info.ord1[0]].imageb[x1_:x2_,y1+j]
                mi = frame[info.ord1[0]].badpixmask[x1_:x2_,y1+j]
                ;si = frame[info.ord1[0]].shift*x1
                ;fi = frame[info.ord1[0]].file
;                FOR k=1,n1-1 do fi=[fi,frame[info.ord1[0]].file]
                for i=1,info.nord1-1 do xi=[xi,xx+frame[info.ord1[i]].shift+shift1]  
                for i=1,info.nord1-1 do yi=[yi,frame[info.ord1[i]].imageb[x1_:x2_,y1+j]]   
                for i=1,info.nord1-1 do mi=[mi,frame[info.ord1[i]].badpixmask[x1_:x2_,y1+j]]   
                for i=1,info.nord1-1 do wi=[wi,var] 
                ;for i=1,info.nord1-1 do si=[si,frame[info.ord1[i]].shift*x1] 
                ;for i=1,info.nord1-1 do FOR k=0,n1-1 do  fi=[fi,frame[info.ord1[i]].file] 
                if info.nord2 gt 0 then begin
                    for i=0,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]    
                    for i=0,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
                    for i=0,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
                    for i=0,info.nord2-1 do wi=[wi,var] 
                 ;   for i=0,info.nord2-1 do si=[si,frame[info.ord2[i]].shift*x1]
                 ;   for i=0,info.nord2-1 do  FOR k=0,n1-1 do fi=[fi,frame[info.ord2[i]].file]
                endif
            ENDELSE
            ;wi = 1./(wi+abs(yi))  ;;; Add positive signal to the variance
            wi = 1./wi  ;; inverse variance
            ind = where(finite(wi) EQ 0, nind) ;; put to zero 
            IF nind GT 0 THEN wi[ind]=0.
            s=sort(xi) & xi=xi[s] & yi=yi[s] & wi=wi[s] & mi=mi[s]; & si=si[s] & fi=fi[s]
            
            ;; Express as flux densities
            yi = yi/info.pixsize ;/dwpix[j]  delta is 1 pixel
            wi = wi * info.pixsize^2
            ;; xi in arcsec relative to the center and distortion corrected
            xi=(xi-pos-dpix[j])*info.pixsize  
            ind = where(abs(xi) LT dw[j], nindw)
            if order eq 1 then spectra.n1[j]=nindw else spectra.n2[j]=nindw
            IF nindw GT 10 THEN BEGIN
                xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind] ; & si=si[ind] & fi=fi[ind]
                ind = where(finite(wi)  AND finite(yi)  AND mi, nind)
                IF nind GT 0 THEN BEGIN
                    xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind]; & si=si[ind] & fi=fi[ind]
 ;                   aa = fltarr(n_elements(xi))+1
 ;                   wtot = [wtot,aa*wpix[j]]
 ;                   stot = [stot,si]
 ;                   ftot = [ftot,fi]
 ;                   xtot = [xtot,xi] 
 ;                   ytot = [ytot,yi]
                    ypsfi = interpol(psf_,xpsf_,xi)
                    den = total(ypsfi^2 * wi)
                    spectrum[j] = total(ypsfi*yi * wi)/den 
                    espectrum[j] = sqrt(1./den)
                    chi2 = (yi-(s>0)*ypsfi)*sqrt(wi) 
                    ;; 5-sigma rejection
                    chi2 = (yi-(spectrum[j]>0)*ypsfi)*sqrt(wi)
                    ind = where((chi2 gt 5. or chi2 LT -3) AND wi GT 0., nind)
                    IF nind GT 0 THEN BEGIN
                        wi[ind]=0.
                        den = total(ypsfi^2 * wi)
                        spectrum[j] = total(ypsfi*yi * wi)/den
                        espectrum[j] = sqrt(1./den)
                    ENDIF
                    ;print,'ratio ',total(yi)/spectrum[j]
                ENDIF ELSE nindw=0
            ENDIF else begin
                spectrum[j]=0
                espectrum[j]=0
            ENDelse
        ENDFOR

;        wtot = wtot[1:*]
;        xtot = xtot[1:*]
;        ytot = ytot[1:*]
;        stot = stot[1:*]
;        ftot = ftot[1:*]
;        dpix_ = dpix[0:y2-y1]
;        save,f='points_'+strtrim(string(order),2)+'.sav',xtot,ytot,wtot,stot,ftot,dpix_,dpix

;; Flux conversion
        fc = calib.cal0[order-1]+calib.cal1[order-1]*(wpix-calib.lref[order-1])+$
             calib.cal2[order-1]*(wpix-calib.lref[order-1])^2
        spectrum = spectrum/(fc*calib.flux_conv[order-1])
        espectrum = espectrum/(fc*calib.flux_conv[order-1])
;; Aperture correction
        IF info.aper THEN BEGIN
            apcorr = interpol(calib.fapcorr,calib.wapcorr,wpix)
            spectrum=spectrum/apcorr
            espectrum=espectrum/apcorr
        ENDIF
    ENDIF ELSE BEGIN
        spectrum[*]=0
        espectrum[*]=0
    ENDELSE

    IF order EQ 1 THEN BEGIN
        spectra.w1 = wpix & spectra.f1a = spectrum & spectra.e1a =espectrum
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        spectra.w2 = wpix & spectra.f2a = spectrum & spectra.e2a =espectrum
    ENDIF
ENDFOR

;; Plot spectra
PlotSpectra

;; Save spectra
SaveSpectra

;; restore limits of profile plot ...
DisplayImage,/profile

END


PRO ExtractSourceVoid, pos
COMMON database

spectra.n = 1
delvarx, points
;;; fit the position
pos=pos*info.sampling

spectra.pos = pos
ComputeCoordinates, pos/info.sampling
spectra.raA = info.ra
spectra.decA = info.dec


;;; Extraction for the two orders
DisplayImage
pos=pos/info.sampling


FOR order=1,2 DO BEGIN

    IF order EQ 1 THEN BEGIN
        shift1=info.shift/info.pixsize & shift2 = 0.
        ima=ima1 
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw1[10:123] 
;            dwpix = info.idw1[10:123] 
        ENDIF  ELSE BEGIN
            wpix = info.iw1[1:126]
;            dwpix = info.idw1[1:126]
        ENDELSE
        spectra.n1[*]=0
        dy_ = interpol(calib.dy1, calib.w1, wpix)
        tot  = total(abs(ima1[pos*info.sampling,*]))
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        shift1=2*info.shift/info.pixsize+info.shift12
        shift2=info.shift/info.pixsize+info.shift12
        ima=ima2
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw2[0:77] 
;            dwpix = info.idw2[0:77] 
        ENDIF ELSE BEGIN
            wpix = info.iw2[0:79]
;            dwpix = info.idw2[0:79]
        ENDELSE
        spectra.n2[*]=0
        dy_ = interpol(calib.dy2, calib.w2, wpix)
        tot = total(abs(ima2[pos*info.sampling,*]),/nan)
    ENDIF
    nx = info.ysize12*info.sampling
    IF info.channel EQ "long" THEN dw=7.172/27.*wpix/2. $
    ELSE  dw=8./12.*wpix/2.
    spectrum = fltarr(n_elements(dw)) 
    espectrum = fltarr(n_elements(dw)) 

    IF tot GT 0 THEN BEGIN
    print,'order tot ',order, tot
        
;;; Plot aperture on the source centered on the profile

        
;;; triple for displaying 
        ypix = findgen(n_elements(dw)*info.sampling)
        dw_= fltarr(n_elements(ypix))
        ii=indgen(n_elements(dw))*info.sampling 
        dw_[ii]=dw & dw_[ii+1]=dw & dw_[ii+2]=dw
        
        IF order EQ 1 THEN wset, info.win_id3 ELSE wset, info.win_id1
    
;;; Compute the flux inside the dw limits by fitting the PSF
;;; to the points ...
        
        indord=where(info.order_mask EQ order, nindord)
        ;print,'no frames in order ',order,' is ', nindord
        x1_ = min(info.x_[indord]) & x2_ = max(info.x_[indord])
        y1  = min(info.y_[indord]) & y2  = max(info.y_[indord])
        ;print,'y1 y2 ',y1,y2
        ypix = findgen(128)+0.5
        IF order EQ 1 THEN BEGIN
            xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
            yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
            IF info.channel EQ "long" THEN ypix = ypix[10:123] ELSE ypix=ypix[1:126]
        ENDIF ELSE IF order EQ 2 THEN BEGIN
            xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
            yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
            IF info.channel EQ "long" THEN ypix = ypix[0:77] ELSE  ypix = ypix[0:79]
        ENDIF
        
        
        s=sort(yc)
        xc=xc[s] & yc=yc[s]
        xpix = interpol(xc,yc,ypix)
        min60 = min(abs(ypix-60.),imin60)
        dpix = xpix - xpix[imin60]
        

        dw = dw * info.pixsize
    
        xx=findgen(x2_-x1_+1)
        x1=fltarr(x2_-x1_+1)+1.
        n1 = n_elements(x1) 

;        ddx = fltarr(n_elements(wpix))*info.magnification*info.sampling
;        wtot=[0]
;        xtot=[0]
;        ytot=[0]
;        stot=[0]
;        ftot=[" "]
        IF info.channel EQ "long" THEN wref=27. ELSE wref=12.
        IF order EQ 1 THEN BEGIN
            wref = psf.c1[0]+psf.c1[1]*wref
            xpsf = psf.x1
            ypsf = psf.y1
        ENDIF ELSE BEGIN
            wref = psf.c2[0]+psf.c2[1]*wref 
            xpsf = psf.x2
            ypsf = psf.y2
        ENDELSE
       
        FOR j=0,y2-y1 DO BEGIN   
            IF j GT 0 THEN BEGIN 
                mag = info.sampling*info.magnification
                plots, col=info.ncolors-1,/dev, lines=0, $
                       ((y2-y1)-[j-1,j])* mag, (pos-dw[j-1:j]/info.pixsize)*mag
                plots, col=info.ncolors-1,/dev, lines=0, $
                       ((y2-y1)-[j-1,j])* mag, (pos+dw[j-1:j]/info.pixsize)*mag
            ENDIF
            IF order EQ 1 THEN BEGIN
                w = psf.c1[0]+psf.c1[1]*wpix[j]
            ENDIF ELSE BEGIN
                w = psf.c2[0]+psf.c2[1]*wpix[j]
            ENDELSE
            xpsf_ = xpsf * w    ;; xpsf_ is in arcsec
            psf_  = ypsf * wref/w  ;; normalized PSF
            psf_ = psf_/dy_[j]  ;; normalize along wavelength direction

            ;; select points
            var = info.rms[x1_:x2_,y1+j]^2
            hot = info.hot[x1_:x2_,y1+j]
            ind=where(hot EQ 1, nind)
            IF nind GT 0 THEN var[ind]=!values.f_nan
            wi = var  ;; inverse variance
            if info.nord1 eq 0 then begin
                xi = xx+frame[info.ord2[0]].shift+shift2 
                yi = frame[info.ord2[0]].imageb[x1_:x2_,y1+j]
                mi = frame[info.ord2[0]].badpixmask[x1_:x2_,y1+j]
                ;si = frame[info.ord2[0]].shift*x1
                ;fi = frame[info.ord2[0]].file
                ;FOR k=1,n1-1 do fi=[fi,frame[info.ord2[0]].file]
                for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]   
                for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
                for i=1,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
                for i=1,info.nord2-1 do wi=[wi,var] 
                ;for i=1,info.nord2-1 do si=[si,frame[info.ord2[i]].shift*x1]
                ;for i=1,info.nord2-1 DO FOR k=0,n1-1 do fi=[fi,frame[info.ord2[i]].file]
            endif else begin
                xi = xx+frame[info.ord1[0]].shift+shift1  
                yi = frame[info.ord1[0]].imageb[x1_:x2_,y1+j]
                mi = frame[info.ord1[0]].badpixmask[x1_:x2_,y1+j]
                ;si = frame[info.ord1[0]].shift*x1
                ;fi = frame[info.ord1[0]].file
;                FOR k=1,n1-1 do fi=[fi,frame[info.ord1[0]].file]
                for i=1,info.nord1-1 do xi=[xi,xx+frame[info.ord1[i]].shift+shift1]  
                for i=1,info.nord1-1 do yi=[yi,frame[info.ord1[i]].imageb[x1_:x2_,y1+j]]   
                for i=1,info.nord1-1 do mi=[mi,frame[info.ord1[i]].badpixmask[x1_:x2_,y1+j]]   
                for i=1,info.nord1-1 do wi=[wi,var] 
                ;for i=1,info.nord1-1 do si=[si,frame[info.ord1[i]].shift*x1] 
                ;for i=1,info.nord1-1 do FOR k=0,n1-1 do  fi=[fi,frame[info.ord1[i]].file] 
                if info.nord2 gt 0 then begin
                    for i=0,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]    
                    for i=0,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
                    for i=0,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
                    for i=0,info.nord2-1 do wi=[wi,var] 
                 ;   for i=0,info.nord2-1 do si=[si,frame[info.ord2[i]].shift*x1]
                 ;   for i=0,info.nord2-1 do  FOR k=0,n1-1 do fi=[fi,frame[info.ord2[i]].file]
                endif
            ENDELSE
            ;wi = 1./(wi+abs(yi))  ;;; Add positive signal to the variance
            wi = 1./wi  ;; inverse variance
            ind = where(finite(wi) EQ 0, nind) ;; put to zero 
            IF nind GT 0 THEN wi[ind]=0.
            s=sort(xi) & xi=xi[s] & yi=yi[s] & wi=wi[s] & mi=mi[s]; & si=si[s] & fi=fi[s]
            
            ;; Express as flux densities
            yi = yi/info.pixsize ;/dwpix[j]  delta is 1 pixel
            wi = wi * info.pixsize^2
            ;; xi in arcsec relative to the center and distortion corrected
            xi=(xi-pos-dpix[j])*info.pixsize  
            ind = where(abs(xi) LT dw[j], nindw)
            if order eq 1 then spectra.n1[j]=nindw else spectra.n2[j]=nindw
            IF nindw GT 10 THEN BEGIN
                xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind] ; & si=si[ind] & fi=fi[ind]
                ind = where(finite(wi)  AND finite(yi)  AND mi, nind)
                IF nind GT 0 THEN BEGIN
                    xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind]; & si=si[ind] & fi=fi[ind]
 ;                   aa = fltarr(n_elements(xi))+1
 ;                   wtot = [wtot,aa*wpix[j]]
 ;                   stot = [stot,si]
 ;                   ftot = [ftot,fi]
 ;                   xtot = [xtot,xi] 
 ;                   ytot = [ytot,yi]
                    ypsfi = interpol(psf_,xpsf_,xi)
                    den = total(ypsfi^2 * wi)
                    spectrum[j] = total(ypsfi*yi * wi)/den 
                    espectrum[j] = sqrt(1./den)
                    chi2 = (yi-(s>0)*ypsfi)*sqrt(wi) 
                    ;; 5-sigma rejection
                    chi2 = (yi-(spectrum[j]>0)*ypsfi)*sqrt(wi)
                    ind = where((chi2 gt 5. or chi2 LT -3) AND wi GT 0., nind)
                    IF nind GT 0 THEN BEGIN
                        wi[ind]=0.
                        den = total(ypsfi^2 * wi)
                        spectrum[j] = total(ypsfi*yi * wi)/den
                        espectrum[j] = sqrt(1./den)
                    ENDIF
                    ;print,'ratio ',total(yi)/spectrum[j]
                ENDIF ELSE nindw=0
            ENDIF else begin
                spectrum[j]=0
                espectrum[j]=0
            ENDelse
        ENDFOR

;        wtot = wtot[1:*]
;        xtot = xtot[1:*]
;        ytot = ytot[1:*]
;        stot = stot[1:*]
;        ftot = ftot[1:*]
;        dpix_ = dpix[0:y2-y1]
;        save,f='points_'+strtrim(string(order),2)+'.sav',xtot,ytot,wtot,stot,ftot,dpix_,dpix

;; Flux conversion
        fc = calib.cal0[order-1]+calib.cal1[order-1]*(wpix-calib.lref[order-1])+$
             calib.cal2[order-1]*(wpix-calib.lref[order-1])^2
        spectrum = spectrum/(fc*calib.flux_conv[order-1])
        espectrum = espectrum/(fc*calib.flux_conv[order-1])
;; Aperture correction
        IF info.aper THEN BEGIN
            apcorr = interpol(calib.fapcorr,calib.wapcorr,wpix)
            spectrum=spectrum/apcorr
            espectrum=espectrum/apcorr
        ENDIF
    ENDIF ELSE BEGIN
        spectrum[*]=0
        espectrum[*]=0
    ENDELSE

    IF order EQ 1 THEN BEGIN
        spectra.w1 = wpix & spectra.f1a = spectrum & spectra.e1a =espectrum
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        spectra.w2 = wpix & spectra.f2a = spectrum & spectra.e2a =espectrum
    ENDIF
ENDFOR

;; Plot spectra
PlotSpectra

;; Save spectra
SaveSpectra

;; restore limits of profile plot ...
DisplayImage,/profile

END





PRO gfunct, X, A, F, pder 
;; Functions with two gaussians

bx1 = EXP(-(X-A[1])^2/A[2]^2)
bx2 = EXP(-(X-A[4])^2/A[5]^2)
F   = A[0]*bx1+A[3]*bx2

;If the procedure is called with four parameters, calculate the 
;partial derivatives. 
IF N_PARAMS() GE 4 THEN $
  pder = [[bx1], [2*A[0]*(X-A[1])* bx1], [2*A[0]*(X-A[1])^2/A[3]^3*bx1],$
          [bx2], [2*A[3]*(X-A[4])* bx2], [2*A[3]*(X-A[4])^2/A[5]^3*bx2]]

END 


PRO Extract2Sources, pos_a, pos_b
COMMON database

spectra.n = 2
delvarx, points

;; order the positions 
IF pos_a GT pos_b THEN BEGIN
    tmp = pos_b
    pos_b = pos_a
    pos_a = tmp
ENDIF

;;; fit the position
pos_a=pos_a*info.sampling
pos_b=pos_b*info.sampling
nx = info.ysize12*info.sampling
profile1=fltarr(nx)
profile2=fltarr(nx)
xprofile = findgen(nx)
IF info.channel EQ "long" THEN BEGIN
    FOR i=0,nx-1 DO profile1[i]=median(ima1[i,30:239])
    FOR i=0,nx-1 DO profile2[i]=median(ima2[i,9:179])
ENDIF ELSE BEGIN
    FOR i=0,nx-1 DO profile1[i]=median(ima1[i,*])
    FOR i=0,nx-1 DO profile2[i]=median(ima2[i,*])   
ENDELSE

;;; Fit a 2 Gaussian profile
;;; check the highest profile
maxi = max([profile1[pos_a]+profile1[pos_b],$
            profile2[pos_a]+profile2[pos_b]],imax)
IF imax EQ 0 THEN profile=profile1 ELSE profile=profile2
pos1=(pos_a-5)>0 & pos2=(pos_b+5)<nx
x = xprofile[pos1:pos2] & y=profile[pos1:pos2]
weights  = replicate(1.,n_elements(x)) 
xmin = min(abs(x-pos_a),ia) & xmin = min(abs(x-pos_b),ib)
a=[y[ia],pos_a,4,y[ib],pos_b,4]
yfit = CURVEFIT(X, Y, weights, A, SIGMA, FUNCTION_NAME='gfunct') 
;;; Oplot the fitted parabola and the limits chosen
wset, info.win_id5
oplot,yfit,x,col=info.ncolors-1
pos_a = a[1] ;;; fitted positions
pos_b = a[4]

;; Compute coordinates
spectra.pos = pos_a
ComputeCoordinates, pos_a/info.sampling
spectra.raA = info.ra
spectra.decA = info.dec
spectra.posB = pos_b
ComputeCoordinates, pos_b/info.sampling
spectra.raB = info.ra
spectra.decB = info.dec


;;; Extraction for the two orders
pos_a=pos_a/info.sampling
pos_b=pos_b/info.sampling
DisplayImage
FOR order=1,2 DO BEGIN
    print,'order,posa ',order,pos_a
    IF order EQ 1 THEN BEGIN
        shift1=info.shift/info.pixsize & shift2 = 0.
        ima=ima1 
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw1[10:123] 
;            dwpix = info.idw1[10:123] 
        ENDIF  ELSE BEGIN
            wpix = info.iw1[1:126]
;            dwpix = info.idw1[1:126]
        ENDELSE
        spectra.n1[*]=0
        dy_ = interpol(calib.dy1, calib.w1, wpix)
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        shift1=2*info.shift/info.pixsize+info.shift12
        shift2=info.shift/info.pixsize+info.shift12
        ima=ima2
        IF info.channel EQ "long" THEN BEGIN
            wpix = info.iw2[0:77] 
;            dwpix = info.idw2[0:77] 
        ENDIF ELSE BEGIN
            wpix = info.iw2[0:79]
;            dwpix = info.idw2[0:79]
        ENDELSE
        spectra.n2[*]=0
        dy_ = interpol(calib.dy2, calib.w2, wpix)
    ENDIF
    nx = info.ysize12*info.sampling
;;; Plot aperture on the source centered on the profile
    IF info.channel EQ "long" THEN dw=7.172/27.*wpix/2. $
    ELSE  dw=8./12.*wpix/2.

;    dw=7.172/27.*wpix/2.  ;; wavelength dependent aperture
;;; triple for displaying 
    ypix = findgen(n_elements(dw)*info.sampling)
    dw_= fltarr(n_elements(ypix))
    ii=indgen(n_elements(dw))*info.sampling
    dw_[ii]=dw & dw_[ii+1]=dw & dw_[ii+2]=dw

    IF order EQ 1 THEN wset, info.win_id3 ELSE wset, info.win_id1
    ;plots, col=info.ncolors-1,reverse(ypix)*info.magnification,$
    ;       (pos_a-dw_)*info.sampling*info.magnification,/dev, lines=0
    ;plots, col=info.ncolors-1,reverse(ypix)*info.magnification,$ 
    ;       (pos_b+dw_)*info.sampling*info.magnification,/dev, lines=0
    
;;; Compute the flux inside the dw limits by fitting the PSF
;;; to the points ...

    indord=where(info.order_mask EQ order, nindord)
    ;print,'no frames in order ',order,' is ', nindord
    x1_ = min(info.x_[indord]) & x2_ = max(info.x_[indord])
    y1  = min(info.y_[indord]) & y2  = max(info.y_[indord])
    ypix = findgen(128)+0.5
    IF order EQ 1 THEN BEGIN
        xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
        yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
        IF info.channel EQ "long" THEN ypix = ypix[10:123]
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
        yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
        IF info.channel EQ "long" THEN ypix = ypix[0:77] ELSE ypix = ypix[0:79]
    ENDIF

    s=sort(yc)
    xc=xc[s] & yc=yc[s]
    xpix = interpol(xc,yc,ypix)
    min60 = min(abs(ypix-60.),imin60)
    dpix = xpix - xpix[imin60]

    IF info.channel EQ "long" THEN dw=7.172/27.*wpix/2.*info.pixsize $
    ELSE  dw=8./12.*wpix/2.*info.pixsize
    xx=findgen(x2_-x1_+1)

    spectrum1 = fltarr(n_elements(dw))  & espectrum1 = fltarr(n_elements(dw)) 
    spectrum2 = fltarr(n_elements(dw))  & espectrum2 = fltarr(n_elements(dw)) 
    FOR j=0,y2-y1 DO BEGIN    

        IF j GT 0 THEN BEGIN 
            mag = info.sampling*info.magnification
            plots, col=info.ncolors-1,/dev, lines=0, $
                   ((y2-y1)-[j-1,j])* mag, (pos_a-dw[j-1:j]/info.pixsize)*mag
            plots, col=info.ncolors-1,/dev, lines=0, $
                   ((y2-y1)-[j-1,j])* mag, (pos_b+dw[j-1:j]/info.pixsize)*mag
        ENDIF


        IF info.channel EQ "long" THEN wref=27. ELSE wref=12.
        IF order EQ 1 THEN BEGIN
            wref = psf.c1[0]+psf.c1[1]*wref
            w = psf.c1[0]+psf.c1[1]*wpix[j]
            xpsf = psf.x1
            ypsf = psf.y1
        ENDIF ELSE BEGIN
            wref = psf.c2[0]+psf.c2[1]*wref 
            w = psf.c2[0]+psf.c2[1]*wpix[j]
            xpsf = psf.x2
            ypsf = psf.y2
        ENDELSE
        xpsf_ = xpsf * w    ;; xpsf_ is in arcsec
        psf_  = ypsf * wref/w  ;; normalized PSF
        psf_ = psf_/dy_[j] ;; normalize along wavelength direction
        ;; select points
        var = info.rms[x1_:x2_,y1+j]^2
        hot = info.hot[x1_:x2_,y1+j]
        ind=where(hot EQ 1, nind)
        IF nind GT 0 THEN var[ind]=!values.f_nan
        wi = var  ;; variance, then transformed in inverse variance
        if info.nord1 eq 0 then begin
            xi = xx+frame[info.ord2[0]].shift+shift2 
            yi = frame[info.ord2[0]].imageb[x1_:x2_,y1+j]
            mi = frame[info.ord2[0]].badpixmask[x1_:x2_,y1+j]
            for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]   
            for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
            for i=1,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
            for i=1,info.nord2-1 do wi=[wi,var] 
        endif else begin
            xi = xx+frame[info.ord1[0]].shift+shift1  
            yi = frame[info.ord1[0]].imageb[x1_:x2_,y1+j]
            mi = frame[info.ord1[0]].badpixmask[x1_:x2_,y1+j]
            for i=1,info.nord1-1 do xi=[xi,xx+frame[info.ord1[i]].shift+shift1]  
            for i=1,info.nord1-1 do yi=[yi,frame[info.ord1[i]].imageb[x1_:x2_,y1+j]]   
            for i=1,info.nord1-1 do mi=[mi,frame[info.ord1[i]].badpixmask[x1_:x2_,y1+j]]   
            for i=1,info.nord1-1 do wi=[wi,var] 
            if info.nord2 gt 0 then begin
                for i=0,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]    
                for i=0,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
                for i=0,info.nord2-1 do mi=[mi,frame[info.ord2[i]].badpixmask[x1_:x2_,y1+j]]   
                for i=0,info.nord2-1 do wi=[wi,var] 
            endif
        ENDELSE
        ;wi = 1./(wi+abs(yi))  ;;; Add positive signal to the variance
        wi = 1./wi  ;; inverse variance
        ind = where(finite(wi) EQ 0, nind) ;; put to zero 
        IF nind GT 0 THEN wi[ind]=0.
        s=sort(xi) & xi=xi[s] & yi=yi[s] & wi=wi[s] & mi=mi[s]

        ;; Express as flux densities
        yi = yi/info.pixsize; /dwpix[j] delta is 1 pixel
        wi = wi * info.pixsize^2

        ;; xi in arcsec distortion corrected
        xi=(xi-dpix[j])*info.pixsize  
        ind=where(xi-pos_a*info.pixsize gt -dw[j] and xi-pos_b*info.pixsize lt dw[j], nindw)

        IF nindw GT 20 THEN BEGIN
            ;print,'nindw ',order,nindw
            xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind]
            ind = where(finite(wi) EQ 1 AND finite(yi) EQ 1 AND mi EQ 1, nind)
            IF nind GT 0 THEN BEGIN
                xi=xi[ind] & yi=yi[ind] & wi=wi[ind] & mi=mi[ind]
                ;; values of psf interpolated
                ypsf1 = interpol(psf_,xpsf_,xi-pos_a*info.pixsize)
                ypsf2 = interpol(psf_,xpsf_,xi-pos_b*info.pixsize)
                ;; fit
                A = total(ypsf1^2*wi)
                B = total(ypsf2^2*wi)
                C = total(ypsf1*ypsf2*wi)
                D = total(ypsf1*yi*wi)
                E = total(ypsf2*yi*wi)
                Den = (A*B-C^2)
                spectrum1[j] = (D*B-C*E)/Den
                spectrum2[j] = (A*E-C*D)/Den
                espectrum1[j] = sqrt(B/Den)
                espectrum2[j] = sqrt(A/Den)
                ;; 5-sigma rejection
                chi2 = (yi-(spectrum1[j]>0)*ypsf1-(spectrum2[j]>0)*ypsf2)*sqrt(wi)
                ind = where((chi2 lt -3. OR chi2 GT 5) AND wi GT 0., nind)
                IF nind GT 0 THEN BEGIN
                    wi[ind]=0.
                    A = total(ypsf1^2*wi)
                    B = total(ypsf2^2*wi)
                    C = total(ypsf1*ypsf2*wi)
                    D = total(ypsf1*yi*wi)
                    E = total(ypsf2*yi*wi)
                    Den = (A*B-C^2)
                    spectrum1[j] = (D*B-C*E)/Den
                    spectrum2[j] = (A*E-C*D)/Den
                    espectrum1[j] = sqrt(B/Den)
                    espectrum2[j] = sqrt(A/Den)
                ENDIF
            ENDIF ELSE nindw=0
        ENDIF ELSE BEGIN
            spectrum1[j]=0 & espectrum1[j]=0
            spectrum2[j]=0 & espectrum2[j]=0
        ENDELSE
    ENDFOR

    ;; Flux conversion
    fc = calib.cal0[order-1]+calib.cal1[order-1]*(wpix-calib.lref[order-1])+$
         calib.cal2[order-1]*(wpix-calib.lref[order-1])^2
    spectrum1 = spectrum1/(fc*calib.flux_conv[order-1])
    espectrum1 = espectrum1/(fc*calib.flux_conv[order-1])
    spectrum2 = spectrum2/(fc*calib.flux_conv[order-1])
    espectrum2 = espectrum2/(fc*calib.flux_conv[order-1])
    ;; Aperture correction
    IF info.aper THEN BEGIN
        apcorr = interpol(calib.fapcorr,calib.wapcorr,wpix)
        spectrum1=spectrum1/apcorr & espectrum1=espectrum1/apcorr
        spectrum2=spectrum2/apcorr & espectrum2=espectrum2/apcorr
    ENDIF

    IF order EQ 1 THEN BEGIN
        spectra.w1 = wpix
        spectra.f1a = spectrum1 & spectra.e1a =espectrum1
        spectra.f1b = spectrum2 & spectra.e1b =espectrum2
    ENDIF ELSE IF order EQ 2 THEN BEGIN
        spectra.w2 = wpix
        spectra.f2a = spectrum1 & spectra.e2a =espectrum1
        spectra.f2b = spectrum2 & spectra.e2b =espectrum2
    ENDIF
ENDFOR

;; Display spectra
PlotSpectra
;; Save spectra
SaveSpectra
;; restore limits of profile plot ...
DisplayImage,/profile


END



PRO MaskImages
COMMON database
;;; Mask sources in the images


;; order 1
shift1=info.shift/info.pixsize & shift2 = 0.
ind = where(info.order_mask EQ 1)
x1_ = min(info.x_[ind]) & x2_ = max(info.x_[ind])
y1  = min(info.y_[ind]) & y2  = max(info.y_[ind])
xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
s=sort(yc)
xc=xc[s] & yc=yc[s]
ypix = findgen(128)
xpix = interpol(xc,yc,ypix)
min60 = min(abs(ypix-60.),imin60)
dpix = xpix - xpix[imin60]
xx=findgen(x2_-x1_+1)
ny = y2-y1+1
; j row, i frame.
ind_ord = where(frame.order EQ 1, nind_ord)

IF info.channel EQ "long" THEN BEGIN
    i1 = 25 & i2 = 56
ENDIF ELSE BEGIN
    i1 = 0 & i2 = 37
ENDELSE 

IF nind_ord GT 0 THEN BEGIN
    limits = info.cliplim[*,0:info.clipn-1]-shift1
    FOR i=0,nind_ord-1 DO BEGIN
        FOR j=0,ny -1 DO BEGIN
            FOR k=0,info.clipn-1 DO begin
                xStart =  x1_+limits[0,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                xEnd   =  x1_+limits[1,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                IF xEnd Gt i1 AND xStart Lt i2 THEN $
                  frame[ind_ord[i]].mask[xStart>i1:xEnd<i2, y1+j] = 0
            ENDFOR
        ENDFOR
    ENDFOR
ENDIF
ind_ord = where(frame.order eq 2, nind_ord)
IF nind_ord GT 0 THEN BEGIN
    limits = info.cliplim[*,0:info.clipn-1]-shift2
    FOR i=0,nind_ord-1 DO BEGIN
        FOR j=0,ny -1 DO BEGIN
            FOR k=0,info.clipn-1 DO begin
                xStart = x1_+limits[0,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                xEnd   = x1_+limits[1,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                IF xEnd Gt i1 AND xStart Lt i2 THEN $
                  frame[ind_ord[i]].mask[xStart>i1:xEnd<i2, y1+j] = 0
        ENDFOR
        ENDFOR
    ENDFOR
ENDIF

;; order 2
shift1=2*info.shift/info.pixsize+info.shift12
shift2=info.shift/info.pixsize+info.shift12
ind = where(info.order_mask EQ 2)
x1_ = min(info.x_[ind]) & x2_ = max(info.x_[ind])
y1  = min(info.y_[ind]) & y2  = max(info.y_[ind])
xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
s=sort(yc)
xc=xc[s] & yc=yc[s]
ypix = findgen(128)
xpix = interpol(xc,yc,ypix)
min60 = min(abs(ypix-60.),imin60)
dpix = xpix - xpix[imin60]
xx=findgen(x2_-x1_+1)
ny = y2-y1+1
ind_ord = where(frame.order EQ 1, nind_ord)
limits = info.cliplim[*,0:info.clipn-1]-shift2
ind_ord = where(frame.order EQ 1, nind_ord)
IF info.channel EQ "long" THEN BEGIN
    i1 = 63 & i2 = 94
ENDIF ELSE BEGIN
    i1 = 42 & i2 =77
ENDELSE

IF nind_ord GT 0 THEN BEGIN
    limits = info.cliplim[*,0:info.clipn-1]-shift1
    FOR i=0,nind_ord-1 DO BEGIN
        FOR j=0,ny -1 DO BEGIN
            FOR k=0,info.clipn-1 DO begin
                xStart =  x1_+limits[0,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                xEnd   =  x1_+limits[1,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                IF xEnd Ge i1 AND xStart Le i2 THEN $
                  frame[ind_ord[i]].mask[xStart:xEnd, y1+j] = 0
            ENDFOR
        ENDFOR
    ENDFOR
ENDIF
ind_ord = where(frame.order eq 2, nind_ord)
IF nind_ord GT 0 THEN BEGIN
    limits = info.cliplim[*,0:info.clipn-1]-shift2
    FOR i=0,nind_ord-1 DO BEGIN
        FOR j=0,ny -1 DO BEGIN
            FOR k=0,info.clipn-1 DO begin
                xStart = x1_+limits[0,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                xEnd   = x1_+limits[1,k]+dpix[y1+j]-frame[ind_ord[i]].shift
                IF xEnd Ge i1 AND xStart Le i2 THEN $
                  frame[ind_ord[i]].mask[xStart:xEnd, y1+j] = 0
        ENDFOR
        ENDFOR
    ENDFOR
ENDIF

END


PRO SubtractBackground
COMMON database

print,'Subtracting background ...'

widget_control,/hourglass
;;; Compute general background and noise

;;; Update masks ...
IF info.clipn GT 0 THEN MaskImages ELSE frame.mask[*,*]=1



;;; take out median image for each order
n = n_elements(frame.order) 
FOR i=0,n-1 DO BEGIN &$
    ind1=where(info.order_mask EQ 1 AND frame[i].mask EQ 1,nind1) &$
    ind2=where(info.order_mask EQ 2 AND frame[i].mask EQ 1,nind2) &$
    ima = frame[i].image &$
    med1 = median(ima(ind1)) &$
    med2 = median(ima(ind2)) &$
    ind1=where(info.order_mask EQ 1,nind1) &$
    ind2=where(info.order_mask EQ 2,nind2) &$
    ima[ind1]=ima[ind1]-med1 &$
    ima[ind2]=ima[ind2]-med2 &$
    frame[i].imageb = ima &$
ENDFOR



;;; Compute background and subtract them for each AOR
aor=frame.aorkey
aor=aor(sort(aor))
aors = aor(uniq(aor))

naor = n_elements(aors)

FOR k=0,naor-1 DO BEGIN &$
    bgr =fltarr(128,128)&$
;   rms =fltarr(128,128)
    indaor = where(frame.aorkey EQ aors[k],nindaor) &$
    for i=0,127 do $
      for j=0,127 do BEGIN &$
        IF info.order_mask[i,j] GT 0 THEN begin &$
          ind = where(frame[indaor].mask[i,j] EQ 1 AND $
                      finite(frame[indaor].imageb[i,j]), nind) &$
          IF nind GT 1 THEN BEGIN &$
;            bbb = mean_ks(frame[indaor[ind]].imageb[i,j],3,/sdev) &$
            bbb = biweight(frame[indaor[ind]].imageb[i,j]) &$
            bgr[i,j] = bbb[0] &$
;           rms[i,j] = bbb[1]/sqrt(nind)
          ENDIF  ELSE bgr[i,j]=!values.f_nan  &$
        ENDIF ELSE bgr[i,j]=median(frame.imageb[i,j])  &$
    ENDFOR &$
;writefits,'rms.fits',rms
    for i=0,nindaor-1 DO frame[indaor[i]].imageb=frame[indaor[i]].imageb-bgr &$
ENDFOR

IF info.out THEN writefits,'bgr.fits',bgr

;;; compute rms from the total ...
rms =fltarr(128,128)
for i=0,127 do $
  for j=0,127 do BEGIN &$
    IF info.order_mask[i,j] GT  0 THEN begin &$
      ind = where(frame.mask[i,j] EQ 1 AND $
                      finite(frame[indaor].imageb[i,j]), nind) &$
      IF nind GT 0 THEN BEGIN &$
         bbb = biweight(frame[indaor[ind]].imageb[i,j]) &$
;        bbb = mean_ks(frame[ind].imageb[i,j],3,/sdev) &$
        rms[i,j] = bbb[1] &$  ;/sqrt(nind)
        ;;; NB: /sqrt(nind) for rms of background, not for single point
      ENDIF  ELSE rms[i,j]=!values.f_nan  &$
  ENDIF  &$
ENDFOR

IF info.out THEN writefits,'rms0.fits',rms

;;; Take out a gradient along the y direction for 1st and 2nd orders
;n = n_elements(frame.order) 
;;;; 1st order
;IF info.channel EQ "long" THEN BEGIN
;    ix1 = 25 & ix2 = 56
;    iy1 = 10 & iy2 = 123
;ENDIF ELSE BEGIN
;    ix1 = 0 & ix2 = 37
;    iy1 = 1 & iy2 = 126
;ENDELSE

;lx1 = (ix2-ix1+1)*(iy2-iy1+1)
;xx = reform(transpose(info.y_[ix1:ix2,iy1:iy2]),lx1)
;x = (findgen(128))[iy1:iy2]

;; don't do any profile subtraction  ?
;
;FOR i=0,n-1 DO BEGIN &$
;    ima=frame[i].imageb &$
;    m1=frame[i].mask & m2=frame[i].badpixmask &$
;    yy = reform(transpose(ima[ix1:ix2,iy1:iy2]),lx1) &$
;    mm1 = reform(transpose(m1[ix1:ix2,iy1:iy2]),lx1) &$
;    mm2 = reform(transpose(m2[ix1:ix2,iy1:iy2]),lx1) &$
;    biw = biweight(yy) &$
;    ind = where(finite(yy) EQ 1 AND abs(yy) LT 5*biw[1] $
;                AND mm1 NE 0 AND mm2 NE 0, nind) &$
;    coeff = poly_fit(xx[ind],yy[ind],3) &$
;    fit = coeff[0]+coeff[1]*x+coeff[2]*x^2+coeff[3]*x^3 &$
;    FOR j=ix1,ix2 DO ima[j,iy1:iy2]=ima[j,iy1:iy2]-fit &$
    ;; subtract gradient
;    frame[i].imageb=ima &$  ;;; don't take out gradient
;ENDFOR

;
;;; 2nd order
;IF info.channel EQ "long" THEN BEGIN
;    ix1 = 63 & ix2 = 94
;    iy1 = 0 & iy2 = 77
;ENDIF ELSE BEGIN
;    ix1 = 42 & ix2 = 77
;    iy1 = 0 & iy2 = 79
;ENDELSE
;lx2=(ix2-ix1+1)*(iy2-iy1+1)
;xx = reform(transpose(info.y_[ix1:ix2,iy1:iy2]),lx2)
;x = (findgen(128))[iy1:iy2]
;ind_ord = where(info.order_mask EQ 0)

;; don't do any profile subtraction  ?
;  or perform something more refined ??

;FOR i=0,n-1 DO BEGIN
;    ima=frame[i].imageb
;    m1=frame[i].mask & m2=frame[i].badpixmask &$
;    yy = reform(transpose(ima[ix1:ix2,iy1:iy2]),lx2)
;    mm1 = reform(transpose(m1[ix1:ix2,iy1:iy2]),lx2) &$
;    mm2 = reform(transpose(m2[ix1:ix2,iy1:iy2]),lx2) &$
;    biw = biweight(yy) &$
;    ind = where(finite(yy) EQ 1 AND abs(yy) LT 5*biw[1] $
;                AND mm1 NE 0 AND mm2 NE 0, nind) &$
;    coeff = poly_fit(xx[ind],yy[ind],3)
;    fit = coeff[0]+coeff[1]*x+coeff[2]*x^2+coeff[3]*x^3
;    FOR j=ix1,ix2 DO ima[j,iy1:iy2]=ima[j,iy1:iy2]-fit
    ;; subtract gradient
;    frame[i].imageb=ima    ;;; don't take away the gradient
;    frame[i].imageb[ind_ord]=!values.f_nan
    ;save,f='prof_'+strtrim(string(i),2)+'.xdr',fit,xx,yy
;ENDFOR


;;; Identify hot pixels
; a) create an image with smoothed rms
rms_smooth = fltarr(128,128)
rms_smooth = rms
rms_smooth(where(info.order_mask EQ 0))=!values.f_nan
rms_smooth = median(rms_smooth,5)
; b) identify pixels
ind = where(rms GT 3*rms_smooth, nind)
info.hot[ind]=1

;;; discard extreme values ... and save a smoothed version of the rms
info.rms = median(rms,3)
IF info.out THEN writefits,'rms.fits',rms


;;; Subtract a running median frame by frame
;print,'subtract a running median ...'
;FOR i=0,n-1 DO BEGIN
;    ima = frame[i].imageb
;    FOR j=0,127 DO BEGIN
;        i1=where(info.order_mask[*,j] EQ 1 AND frame[i].mask[*,j] $
;                 AND info.hot[*,j] EQ 0, n1)
;        IF n1 GT 0 THEN BEGIN
;            ind1=where(info.order_mask[*,j] EQ 1)
;            ima[ind1,j]=ima[ind1,j]-median(ima[i1,j])
;        ENDIF
;        i2=where(info.order_mask[*,j] EQ 2 AND frame[i].mask[*,j] $
;                 AND info.hot[*,j] EQ 0, n2)
;        IF n2 GT 0 THEN BEGIN
;            ind2=where(info.order_mask[*,j] EQ 2)
;            ima[ind2,j]=ima[ind2,j]-median(ima[i2,j])
;        ENDIF
;    ENDFOR
;    frame[i].imageb = ima
;ENDFOR



;;; Alternatively, we can group frames with the same shift and order 
;;; and compute a running median. We should have a better statistics
;;; in this way.

print, 'subtract a running median ...'
FOR order=1,2 DO BEGIN
    indo = where(frame.order EQ order, nindo)
    IF nindo GT 0 THEN BEGIN
        shift = frame[indo].shift
        s = round(shift[sort(shift)])
        s = s[uniq(s)]
        FOR l=0,n_elements(s)-1 DO BEGIN
            ind = where(frame.order EQ order AND round(frame.shift) EQ s[l], nind)
            FOR j=0,127 DO BEGIN
                p=[0]
                i1=where(info.order_mask[*,j] EQ 1 AND frame[ind[0]].mask[*,j] $
                         AND info.hot[*,j] EQ 0, n1)
                FOR k=0,nind-1 DO IF n1 GT 0 THEN p=[p,frame[ind[k]].imageb[i1,j]]
                IF n_elements(p) GT 2 THEN BEGIN
                    p=p[1:*]
                    medp=median(p)
                    ind1=where(info.order_mask[*,j] EQ 1)
                    FOR k=0,nind-1 DO $
                      frame[ind[k]].imageb[ind1,j]=frame[ind[k]].imageb[ind1,j]-medp
                ENDIF
                p=[0]
                i1=where(info.order_mask[*,j] EQ 2 AND frame[ind[0]].mask[*,j] $
                         AND info.hot[*,j] EQ 0, n1)
                FOR k=0,nind-1 DO IF n1 GT 0 THEN p=[p,frame[ind[k]].imageb[i1,j]]
                IF n_elements(p) GT 2 THEN BEGIN
                    p=p[1:*]
                    medp=median(p)
                    ind1=where(info.order_mask[*,j] EQ 2)
                    FOR k=0,nind-1 DO $
                      frame[ind[k]].imageb[ind1,j]=frame[ind[k]].imageb[ind1,j]-medp
                ENDIF
            ENDFOR
        ENDFOR
    ENDIF
ENDFOR





;;; Save images background subtracted
frameName = ['a','b','c','d','e','f','g','h','i']
FOR order=1,2 DO BEGIN
    indo = where(frame.order EQ order, nindo)
    IF nindo GT 0 THEN BEGIN
        shift = frame[indo].shift
        s = round(shift[sort(shift)])
        s = s[uniq(s)]
        FOR l=0,n_elements(s)-1 DO BEGIN
            ind = where(frame.order EQ order AND round(frame.shift) EQ s[l], nind)
            image = fltarr(128,128)
            FOR i=0,127 DO FOR j=0,127 DO BEGIN
                m = frame[ind].badpixmask[i,j] & im = frame[ind].imageb[i,j]
                idm = where(m EQ 1, nidm)
                IF nidm GT 0 THEN image[i,j]= mean(im[idm],/nan) $
                ELSE image[i,j]=!values.f_nan
;                IF nidm GT 0 AND finite(image[i,j]) EQ 0 THEN print,'nan ',i,j
            ENDFOR 
            ;writefits,'image.fits',image
;            IF order EQ 1 THEN stop
            h = headfits(frame[ind[0]].file)
            ;; interpolate along column hot pixels ...
            indhot = where(info.hot EQ 1 or image LT -3*rms, nindhot)
            IF nindhot GT 0 THEN $
              FOR i=0,nindhot-1 DO BEGIN &$
                ix = indhot[i] mod 128 & iy = indhot[i] / 128 &$
                IF iy EQ 0  OR iy EQ 127 THEN $
                  image[ix,iy] = !values.f_nan ELSE BEGIN &$
                    IF info.hot[ix,iy-1] EQ 1 AND info.hot[ix,iy+1] EQ 1 THEN $
                      image[ix,iy] = !values.f_nan ELSE $
                        image[ix,iy] = (image[ix,iy-1]+image[ix,iy+1])/2. &$
                ENDELSE &$
            ENDFOR
            ind = where(info.order_mask EQ 0)
            image[ind]=!values.f_nan
            IF info.out THEN BEGIN
                print,'Saving images ..'
                ord = strtrim(string(order),2)
                writefits,'ima_'+ord+frameName[l]+'.fits',image, h &$
                ind = where(finite(image) EQ 0,nind) 
                bmask = info.bmask
                IF nind GT 0 THEN bmask[ind]=bmask[ind]+2^11
                writefits,'mask_'+ord+frameName[l]+'.fits',bmask, h 
                ind = where(frame.order EQ order AND round(frame.shift) EQ s[l], nind)                
                rms_image = sqrt(info.rms^2+abs(image))/sqrt(nind)
                ind = where(finite(rms_image) EQ 0, nind)
                IF nind GT 0 THEN rms_image[ind] = max(rms_image,/nan)
                writefits,'rms_'+ord+frameName[l]+'.fits',rms_image, h 
            ENDIF
        ENDFOR
    ENDIF
ENDFOR

;;; Interpolate frames to have the same wavelength as the
;;; central ridge (skipping hot pixels ...)

;;; order 1
;a) find out the wavelengths of central ridge pixels
xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
w  = calib.w1
s=sort(yc) & xc=xc[s] & yc=yc[s] & w=w[s]

IF info.channel EQ "long" THEN BEGIN
    ix1 = 25 & ix2 = 56
    iy1 = 10 & iy2 = 123
ENDIF ELSE BEGIN
    ix1 = 0 & ix2 = 37
    iy1 = 1 & iy2 = 126
ENDELSE
ypix = iy1+findgen(iy2-iy1+1)+0.5
xpix = interpol(xc,yc,ypix)
wpix = interpol(w,yc,ypix)
info.iw1[ypix-0.5]=wpix

;b) substitute in imageb fluxes corresponding to the same wavelengths
;   of the central ridge in the pixels of same rows.
FOR i=0,n_elements(frame.shift)-1 DO BEGIN
    FOR j=-19,19 DO BEGIN
        xx = xpix+j
        ind = where(xx GE ix1 AND xx LE ix2, nind)
        IF nind GT 0 THEN BEGIN
            xx_ = xx[ind] & ypix_ =ypix[ind]
            ww = calib.wave[xx_,ypix_]
            image = frame[i].imageb[xx_,ypix_]
            hot = info.hot[xx_,ypix_]
            indh=where(hot EQ 0 AND finite(image) EQ 1, nindh)
            image_interp = fltarr(n_elements(image)) 
            IF nindh GT 1 THEN $
              image_interp = interpol(image[indh],ww[indh],wpix[ind])
            indh=where(hot EQ 1 or finite(image) EQ 0, nindh)
            IF nindh GT 0 THEN image_interp[indh]=!values.f_nan
            frame[i].imageb[xx_,ypix_]=image_interp
        ENDIF
    endfor
ENDFOR

;;; order 2
;a) find out the wavelengths of central ridge pixels
xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
w  = calib.w2
s=sort(yc) & xc=xc[s] & yc=yc[s] & w=w[s]

IF info.channel EQ "long" THEN BEGIN
    ix1 = 63 & ix2 = 94
    iy1 = 0 & iy2 = 77
ENDIF ELSE BEGIN
    ix1 = 42 & ix2 = 77
    iy1 = 0 & iy2 = 79
ENDELSE
ypix = iy1+0.5+findgen(iy2-iy1+1)
xpix = interpol(xc,yc,ypix)
wpix = interpol(w,yc,ypix)
info.iw2[ypix-0.5]=wpix

;b) substitute in imageb fluxes corresponding to the same wavelengths
;   of the central ridge in the pixels of same rows.
FOR i=0,n_elements(frame.shift)-1 DO BEGIN
    FOR j=-19,19 DO BEGIN
        xx = xpix+j
        ind = where(xx GE ix1 AND xx LE ix2, nind)
        IF nind GT 0 THEN BEGIN
            xx_ = xx[ind] & ypix_ =ypix[ind]
            ww = calib.wave[xx_,ypix_]
            image = frame[i].imageb[xx_,ypix_]
            hot = info.hot[xx_,ypix_]
            indh=where(hot EQ 0 AND finite(image) EQ 1, nindh)
            image_interp = fltarr(n_elements(image)) 
            IF nindh GT 1 THEN $
              image_interp = interpol(image[indh],ww[indh],wpix[ind])
            indh=where(hot EQ 1 or finite(image) EQ 0, nindh)
            IF nindh GT 0 THEN image_interp[indh]=!values.f_nan
            frame[i].imageb[xx_,ypix_]=image_interp
        ENDIF
    endfor
ENDFOR

;print,' Saving the frame '
IF info.out THEN save,f='cube.xdr',frame,info



END



PRO ComputeOrder, order
COMMON database

print,'Computing order ',order

;;; compute xp1,yp1,w1 in the center of the order 
IF order EQ 1 THEN BEGIN
    xc = 0.5 * (calib.x_start[calib.ind1]+calib.x_end[calib.ind1])
    yc = 0.5 * (calib.y_start[calib.ind1]+calib.y_end[calib.ind1])
    w  = calib.w1
ENDIF ELSE IF order EQ 2 THEN BEGIN
    xc = 0.5 * (calib.x_start[calib.ind2]+calib.x_end[calib.ind2])
    yc = 0.5 * (calib.y_start[calib.ind2]+calib.y_end[calib.ind2])
    w  = calib.w2
ENDIF
    

s=sort(yc)
xc=xc[s] & yc=yc[s] & w=w[s]
ypix = findgen(128)+0.5
xpix = interpol(xc,yc,ypix)
;;; minimum distortion around ypix = 60 for 1st and 2nd order LL
min60 = min(abs(ypix-60.),imin60)
dpix = xpix - xpix[imin60]
wpix = interpol(w,yc,ypix)

ind=where(info.order_mask EQ order)
x1_ = min(info.x_[ind]) & x2_ = max(info.x_[ind])
y1  = min(info.y_[ind]) & y2  = max(info.y_[ind])

;;; Factor 3 only for visualization ... it is not better to include it
;;; in info.magnification ? No, the factor 3 is an oversampling factor
;;; in the X coordinate, now info.sampling
xx=findgen(x2_-x1_+1)
nx = (round(max(frame.shift))+max(xx)+info.shift/info.pixsize)*info.sampling

xdata = findgen(nx)/info.sampling
ny = y2-y1+1

IF order EQ 1 THEN  $
  outima = fltarr(info.ysize12*info.sampling,info.xsize2*info.sampling) $
ELSE $
  outima = fltarr(info.ysize12*info.sampling,info.xsize1*info.sampling)

;;;; CHECK AGAIN ....
IF order EQ 1 THEN BEGIN
    shift1=info.shift/info.pixsize & shift2 = 0.
ENDIF ELSE if order EQ 2 THEN BEGIN
    shift1=2*info.shift/info.pixsize+info.shift12
    shift2=info.shift/info.pixsize+info.shift12
endif

nx = n_elements(xx) 
for j=0,y2-y1 do begin    
    if info.nord1 eq 0 then begin
        xi = xx+frame[info.ord2[0]].shift+shift2 
        yi = frame[info.ord2[0]].imageb[x1_:x2_,y1+j]
        mi = xx+x1_ 
        mf = replicate(info.ord2[0],nx)
        for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]   
        for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
        for i=1,info.nord2-1 do mi=[mi,xx+x1_]   
        for i=1,info.nord2-1 do mf=[mf,replicate(info.ord2[i],nx)]   
    endif else begin
        xi = xx+frame[info.ord1[0]].shift+shift1  
        yi = frame[info.ord1[0]].imageb[x1_:x2_,y1+j]
        mi = xx+x1_ 
        mf = replicate(info.ord1[0],nx)
        for i=1,info.nord1-1 do xi=[xi,xx+frame[info.ord1[i]].shift+shift1]  
        for i=1,info.nord1-1 do yi=[yi,frame[info.ord1[i]].imageb[x1_:x2_,y1+j]]   
        for i=1,info.nord1-1 do mi=[mi,xx+x1_]   
        for i=1,info.nord1-1 do mf=[mf,replicate(info.ord1[i],nx)]   
        if info.nord2 gt 0 then begin
            for i=1,info.nord2-1 do xi=[xi,xx+frame[info.ord2[i]].shift+shift2]    
            for i=1,info.nord2-1 do yi=[yi,frame[info.ord2[i]].imageb[x1_:x2_,y1+j]]   
            for i=1,info.nord2-1 do mi=[mi,xx+x1_]   
            for i=1,info.nord2-1 do mf=[mf,replicate(info.ord2[i],nx)]   
        endif
    endelse
    ;IF order EQ 1 AND j EQ 2 THEN stop

;;; smoothing with a b-cubic spline   
;; wavelength (pixel position)
    y=xi-dpix[j]
    s=sort(y)  & xi=xi[s] & yi=yi[s] & mi=mi[s] & mf=mf[s]  
    y=y[s] & f=yi   
    ;; knots
    f0=f
    ind = where(finite(f))
    x_ = y[ind]
    dx = x_[1:*]-x_[0:n_elements(x_)-2] 
    indgap = where(dx GE 0.5,ngap)
    xg1 = x_[indgap+1] & xg0 = x_[indgap]
    ;gap = median(xg1-xg0)
    dx = xg0[1:*]-xg1[0:n_elements(xg1)-2]
    fullbkpt = [min(x_),xg1]
    ind = where(dx GT 2., nind)
    IF nind GT 0 THEN $
      FOR i=0,nind-1 DO $
        fullbkpt = [fullbkpt,xg1[ind[i]]+1.+findgen(floor(dx[ind[i]]))]
    fullbkpt=fullbkpt[sort(fullbkpt)]


;; weight
    ind = where(finite(f0))  
    biw = biweight(f[ind])   
    invvar = fltarr(n_elements(f))+1./biw[1]   
    ;invvar = fltarr(n_elements(f))+1.
; deleting extreme outliers
    ;ind = where(abs(f-median(f)) gt 20 * biw[1], nind)   
    ;if nind gt 0 then invvar[ind]=0.   
; deleting nans
    ind = where(finite(f) eq 0, nind)   
    if nind gt 0 then invvar[ind]=0.   
    if nind gt 0 then f[ind]=0.   
;; run the code
    sset = 0   
    error_code = bspline_fit(y,f, invvar, sset, fullbkpt=fullbkpt, yfit=yfit)   
;plot,y,f,ps=3,/xs
;oplot,y,yfit
;oplot,fullbkpt,bspline_valu(fullbkpt,sset),ps=4

;; residuals
    ind = where(invvar gt 0)
    res = f-yfit   
    biw = biweight(res[ind])
    disp = biw[1]+sqrt(yfit>0)
    indout = where(res gt 7*disp OR res LT -7*biw[1], nindout)   
    if nindout gt 0 then invvar[indout]=0   
    ;;; check frame.mask
    if nindout gt 0 then begin &$
        mmi = mi [indout] &$
        mmf = mf [indout] &$
        ind = where(mmi gt 0,nind) &$
        if nind gt 0 then begin &$
            for l=0,nind-1 do frame[mmf[l]].badpixmask[mmi[l],y1+j]=0 &$
        endif &$
    ENDIF

    ind = where(invvar GT 0)
    x_ = y[ind]
    dx = x_[1:*]-x_[0:n_elements(x_)-2] 
    indgap = where(dx GE 0.5,ngap)
    xg1 = x_[indgap+1] & xg0 = x_[indgap]
    ;gap = median(xg1-xg0)
    dx = xg0[1:*]-xg1[0:n_elements(xg1)-2]
    fullbkpt = [min(x_),xg1]
    ind = where(dx GT 2., nind)
    IF nind GT 0 THEN $
      FOR i=0,nind-1 DO $
        fullbkpt = [fullbkpt,xg1[ind[i]]+1.+findgen(floor(dx[ind[i]]))]
    fullbkpt=fullbkpt[sort(fullbkpt)]

    sset = 0   
    error_code = bspline_fit(y,f, invvar, sset, fullbkpt=fullbkpt, yfit=yfit)   
;; second tour
    ind = where(invvar gt 0)   
    res = f-yfit   
    biw = biweight(res[ind])   
    disp = biw[1]+sqrt(yfit>0)
    indout = where(res gt 4*disp OR res LT -3*biw[1], nindout)   
    if nindout gt 0 then invvar[indout]=0   
    ;;; check frame.mask
    ;if nindout gt 0 then begin
    ;    mmi = mi [indout]
    ;    mmf = mf [indout] &$
    ;    ind = where(mmi gt 0,nind)
    ;    if nind gt 0 then begin
    ;        for l=0,nind-1 do frame[mmf[l]].badpixmask[mmi[l],y1+j]=0
    ;    endif
    ;endif
    ind = where(invvar GT 0)
    x_ = y[ind]
    dx = x_[1:*]-x_[0:n_elements(x_)-2] 
    indgap = where(dx GE 0.5,ngap)
    xg1 = x_[indgap+1] & xg0 = x_[indgap]
    ;gap = median(xg1-xg0)
    dx = xg0[1:*]-xg1[0:n_elements(xg1)-2]
    fullbkpt = [min(x_),xg1]
    ind = where(dx GT 2., nind)
    IF nind GT 0 THEN $
      FOR i=0,nind-1 DO $
        fullbkpt = [fullbkpt,xg1[ind[i]]+1.+findgen(floor(dx[ind[i]]))]
    fullbkpt=fullbkpt[sort(fullbkpt)]

    sset = 0   
    error_code = bspline_fit(y,f, invvar, sset, fullbkpt=fullbkpt, yfit=yfit)   
    IF error_code EQ -1 THEN print,order,j,error_code
;a) compute spine values from coeffs
    ind = where(finite(f0))    
    n = max(y[ind])-min(y[ind])
    xdata = round(min(y))+findgen(n*info.sampling)/info.sampling
    xdata_ = xdata;+dpix[y1+j];;+0.5/info.sampling ;;; add half pixel ?
;    IF order EQ 2 THEN xdata_ = xdata_ - info.shift12
    ind = where(xdata_ GT min(y) AND xdata_ LT max(y), nind)
    xdata=xdata[ind] & xdata_=xdata_[ind]
    yfit = bspline_valu(xdata_, sset) 
    xp = round(min(y))*info.sampling+findgen(n*info.sampling)
    xp = xp[ind]
    FOR k=0,info.sampling-1 DO $
      outima[xp,info.sampling*j+k]=yfit   
endfor

if order eq 1 then ima1 = outima else ima2=outima
;save,file='images.xdr',ima1,ima2

end


pro DisplayImage, profile=profile
common database


ima1_ =  reverse(transpose(ima1)) ;; exchange x and y axis
ima2_ =  reverse(transpose(ima2)) ;; exchange X and Y axis

IF NOT keyword_set(profile) THEN BEGIN
;;; order 1
    wset, info.win_id3
    xsize = info.magnification*info.xsize2*info.sampling
    ysize = info.magnification*info.ysize12*info.sampling
    ima_ = congrid(ima1_,xsize,ysize)<info.maxscale>info.minscale
    ima_ = (info.ncolors-2) * (ima_ - info.minscale) / $
           (info.maxscale-info.minscale)
    tv,ima_
    
;;; ORDER 2
    ima2_ =  reverse(transpose(ima2)) ;; exchange X and Y axis
    wset, info.win_id1
    xsize = info.magnification*info.xsize1*info.sampling
    ysize = info.magnification*info.ysize12*info.sampling
    ima_ = congrid(ima2_,xsize,ysize)<info.maxscale>info.minscale
    ima_ = (info.ncolors-2) * (ima_ - info.minscale) / $
           (info.maxscale-info.minscale)
    tv,ima_
ENDIF


;; Plot the profiles

ny = info.ysize12*info.sampling
profile1=fltarr(ny)
profile2=fltarr(ny)
yprofile = findgen(ny)

FOR i=0,ny-1 DO profile1[i]=median(ima1_[30:239,i])
FOR i=0,ny-1 DO profile2[i]=median(ima2_[9:179,i])

wset, info.win_id5
;print,'x limits are ',info.maxscale,info.minscale
plot,xma=[0,0],yma=[0,0],[0,0],[0,ny],/xs,/ys,yr=[0,ny],xr=[info.maxscale,info.minscale],col=info.ncolors-2


IF info.clipn GT 0 THEN $
  FOR i=0,info.clipn-1 DO BEGIN
    x1 = info.maxscale & x2=info.minscale
    l1 = info.cliplim[0,i]*info.sampling & l2 = info.cliplim[1,i]*info.sampling
    polyfill, [x1,x2,x2,x1,x1], [l1,l1,l2,l2,l1],color=info.ncolors/5,noclip=0
ENDFOR

oplot,profile1,yprofile,col=info.ncolors-1
oplot,profile2,yprofile,col=info.ncolors-2,lines=2

;; Overplot the center ...
oplot,[info.maxscale,info.minscale], info.x0*info.sampling*[1,1]
;save,f='test.save',ima1_,ima2_

END




PRO IrsLow_event, ev
COMMON database

widget_control, /draw_button_events, get_uvalue=type, ev.id
thisEvent = Tag_Names(ev, /Structure_Name)

IF thisEvent EQ "WIDGET_TRACKING" THEN $
  widget_control, ev.id, /INPUT_FOCUS

IF thisEvent NE "WIDGET_TRACKING" THEN BEGIN
    IF type EQ "SCALE_MIN" THEN BEGIN
        widget_control, info.text[0], get_value=tmp
        datatype=strnumber(tmp[0],val)
        IF datatype THEN BEGIN
            IF val lt info.maxscale THEN BEGIN
                info.minscale=val & DisplayImage 
            ENDIF ELSE $
                  widget_control, info.text[0], set_value=strtrim(info.minscale,2)
        ENDIF ELSE $
          widget_control, info.text[0], set_value=strtrim(info.maxscale,2)
    ENDIF ELSE IF type EQ "SCALE_MAX" THEN BEGIN
        widget_control, info.text[1], get_value=tmp
        datatype=strnumber(tmp[0],val)
        IF datatype THEN BEGIN
            IF val gt info.minscale THEN BEGIN
                info.maxscale=val & DisplayImage
            ENDIF ELSE $
              widget_control, info.text[1], set_value=strtrim(info.maxscale,2)
        ENDIF ELSE $
          widget_control, info.text[1], set_value=strtrim(info.minscale,2)
    ENDIF ELSE IF type EQ "DRAW_WINDOW1" OR type EQ "DRAW_WINDOW3" OR  $
      type EQ "DRAW_WINDOW5" AND ev.press THEN BEGIN
        ;; Get cursor, compute coordinate
        cursor=[ev.x, ev.y]/info.sampling/info.magnification
        print,'cursor ',cursor[1]
        ComputeCoordinates,cursor[1]
        ;; Keyboard options
        opt=""
        idx = Where(Byte('qQ') EQ ev.ch, cnt) & IF cnt THEN opt = "QUIT"
        idx = Where(Byte('bB') EQ ev.ch, cnt) & IF cnt THEN opt = "BACKGROUND"
        idx = Where(Byte('a') EQ ev.ch, cnt) & IF cnt THEN opt = "ALIGN"
        idx = Where(Byte('c') EQ ev.ch, cnt) & IF cnt THEN opt = "CLIP"
        idx = Where(Byte('C') EQ ev.ch, cnt) & IF cnt THEN opt = "UNCLIP"
        idx = Where(Byte('x') EQ ev.ch, cnt) & IF cnt THEN opt = "EXTRACT"
        idx = Where(Byte('X') EQ ev.ch, cnt) & IF cnt THEN opt = "EXTRACT"
        idx = Where(Byte('2') EQ ev.ch, cnt)  & IF cnt THEN opt = "EXTRACT2"
        idx = Where(Byte('?hH') EQ ev.ch, cnt) & IF cnt THEN opt = "HELP"
        CASE opt OF 
            "HELP": xdisplayfile, 'dummy',title='Key actions', $
              group=ev.top, height=15, width=80, $ ;;Help
              text=['?hH: this help', $
                    'a:  click a on 1st and 2nd order sources to align their profiles',$
                    'c:  click c on left and right limits of the profile to clip',$
                    'C:  click C on left and right limits of the clipped profile to unclip',$
                    'bB:  recompute the background (do only if more than 2 observations)',$
                    'xX:  click x on the profile to extract a spectrum',$
                    '2:   click on two nearby sources to extract two spectra',$
                    '  click on the spectrum to see points and PSF fitted',$
                    'qQ:  quit']
            "EXTRACT": BEGIN
                wind=0
                IF type EQ "DRAW_WINDOW1" THEN wind=2 ELSE IF $
                  type EQ "DRAW_WINDOW3" then wind=1
                ExtractSource,cursor[1],wind
            END
            "EXTRACTVOID": ExtractSourceVoid,cursor[1]
            "EXTRACT2": BEGIN ;; mark a source, display on window5
                print,'info.action ',info.action
                IF info.action EQ "NONE" THEN BEGIN
                    info.action = "EX2"
                    info.tmp = cursor[1]
                    print,'cursor is ',cursor[0]
                ENDIF ELSE IF info.action EQ "EX2" THEN BEGIN
                    info.action = "NONE"
                    print,'sources at: ',info.tmp, cursor[1]
                    Extract2Sources, info.tmp, cursor[1]
                ENDIF
            END
            "ALIGN": BEGIN
                IF info.action EQ "NONE" THEN BEGIN
                    info.action = "ALIGN"
                    ;; save position
                    info.tmp = cursor[1]
                ENDIF ELSE IF info.action EQ "ALIGN" THEN BEGIN
                    info.action = "NONE"
                    ;; save position and detect order                    
                    IF type EQ "DRAW_WINDOW1" THEN BEGIN
                        pos1 = info.tmp
                        pos2 = cursor[1]
                    ENDIF ELSE IF type EQ "DRAW_WINDOW3" THEN BEGIN
                        pos1 = cursor[1]
                        pos2 = info.tmp
                    ENDIF ELSE BEGIN
                        print,'Please click on the image windows !'
                    ENDELSE
                    AlignOrders, pos1, pos2
                ENDIF
            END
            "CLIP": BEGIN 
                IF info.action EQ "NONE" THEN BEGIN
                    info.action = "CLIP"
                    info.clipn =info.clipn+1 
                    info.cliplim[0,info.clipn-1] = cursor[1]
                ENDIF ELSE IF info.action EQ "CLIP" THEN BEGIN
                    info.action = "NONE"
                    info.cliplim[1,info.clipn-1] = cursor[1]
                    ;; Order clip limits
                    IF info.cliplim[0,info.clipn-1] GT $
                      info.cliplim[1,info.clipn-1]  THEN BEGIN
                        tmp = info.cliplim[1,info.clipn-1]
                        info.cliplim[1,info.clipn-1]=info.cliplim[0,info.clipn-1]
                        info.cliplim[0,info.clipn-1]=tmp
                    ENDIF
                    ;;; check limits
                    ny = info.ysize12*info.sampling
                    profile=fltarr(ny)
                    yprofile = findgen(ny)
                    FOR i=0,ny-1 DO $
                      profile[i]=median(ima1[i,30:239])+median(ima2[i,9:179])
                    ind=where(profile GT 0, nind)
                    imin = min(ind)/info.sampling & imax = max(ind)/info.sampling
                    IF info.cliplim[0,info.clipn-1] LT imin THEN $
                      info.cliplim[0,info.clipn-1] = imin
                    IF info.cliplim[1,info.clipn-1] GT imax THEN $
                      info.cliplim[1,info.clipn-1] = imax
                    print, imin, imax
                    print, info.cliplim[*,info.clipn-1]
                    DisplayImage
                    wind=[info.win_id1,info.win_id3]
                    FOR k=0,1 DO BEGIN
                        wset, wind[k]
                        FOR i=0,info.clipn-1 DO BEGIN
                            plots, col=info.ncolors-1, $
                                   [0,info.xsize2*info.sampling*info.magnification],$
                                   info.cliplim[0,i]*[1,1] $
                                   *info.sampling*info.magnification,/dev, lines=2
                            plots, col=info.ncolors-1,  $
                                   [0,info.xsize2*info.sampling*info.magnification],$
                                   info.cliplim[1,i]*[1,1] $
                                   *info.sampling*info.magnification,/dev, lines=2
                        ENDFOR
                    ENDFOR
                ENDIF 
            END
            "UNCLIP": BEGIN 
                IF info.action EQ "NONE" THEN BEGIN
                    info.action = "UNCLIP"
                    info.unclip = cursor[1]
                ENDIF ELSE IF info.action EQ "UNCLIP" THEN BEGIN
                    info.action = "NONE"
                    unclip1=info.unclip & unclip2 = cursor[1]
                    info.unclip = 0.
                    ;; order unclip limits
                    IF unclip2 LT unclip1 THEN BEGIN
                        tmp = unclip2
                        unclip2=unclip1
                        unclip1=tmp
                    ENDIF
                    ;; search for clipping limits
                    IF info.clipn GT 0 THEN BEGIN
                        ind = where(info.cliplim[0,0:info.clipn-1] GE unclip1 $
                                    AND info.cliplim[1,0:info.clipn-1] LE $
                                    unclip2, nind)
                        IF nind GT 0 THEN BEGIN
                            info.clipn = info.clipn-nind
                            info.cliplim[0,ind]=0.
                            info.cliplim[1,ind]=0.
                            s=reverse(sort(info.cliplim[0,*]))
                            info.cliplim[0,*]=info.cliplim[0,s]
                            info.cliplim[1,*]=info.cliplim[1,s]
                        ENDIF
                    ENDIF
                    DisplayImage
                ENDIF
            END
            "BACKGROUND": BEGIN
                SubtractBackground
                ComputeOrder, 1
                ComputeOrder, 2
                DisplayImage
            END
            "QUIT": BEGIN
                widget_control, /reset
            END
            ELSE:
        ENDCASE
    ENDIF ELSE IF type EQ "DRAW_WINDOW2"  AND ev.press $
      AND spectra.n GT 0 THEN BEGIN
    ;;; Events in the spectral windows 
    ;;; 2nd order
        ;; Get x position
        j=ev.x/info.sampling/info.magnification+0.5
        j=round(j)
        j=(info.xsize1-j)<(info.xsize1-1)>0
        ;print,'j ',ev.x,j
        SelectPoints, j, 2
    ENDIF ELSE IF type EQ "DRAW_WINDOW4"  AND ev.press $
      AND spectra.n GT 0 THEN BEGIN
    ;;; 1st order
        ;; Get x position
        j=ev.x/info.sampling/info.magnification+0.5
        ;print,' ',j,round(j)
        j=round(j)
        j=(info.xsize2-j)<(info.xsize2-1)>0
        ;print,'j ',info.sampling*info.magnification,ev.x,j
        SelectPoints, j, 1
    ENDIF ELSE IF type EQ "DRAW_WINDOW6" AND ev.press THEN BEGIN 
    ;;; Events in the fitting window
        ;; Get cursor
        cursor=(convert_coord([ev.x, ev.y],/dev))[0:1]
        print,'cursor ',cursor
        ;; Keyboard options
        opt=""
        idx = Where(Byte('?hH') EQ ev.ch, cnt) & IF cnt THEN opt = "HELP"
        idx = Where(Byte('dD') EQ ev.ch, cnt) & IF cnt THEN opt = "DELETE"
        idx = Where(Byte('uU') EQ ev.ch, cnt) & IF cnt THEN opt = "UNDELETE"
        idx = Where(Byte('fF') EQ ev.ch, cnt) & IF cnt THEN opt = "FIT"
        idx = Where(Byte('sS') EQ ev.ch, cnt) & IF cnt THEN opt = "SAVE"
        idx = Where(Byte('w') EQ ev.ch, cnt) & IF cnt THEN opt = "RESIZE"
        idx = Where(Byte('W') EQ ev.ch, cnt) & IF cnt THEN opt = "ORIGINALSIZE"
        CASE opt OF
            "HELP": xdisplayfile, 'dummy',title='Key actions', $
              group=ev.top, height=15, width=80, $ ;;Help
              text=['?hH: this help', $
                    'dD: delete a point', $
                    'uU: undelete a point', $
                    'fF: redo the fit', $
                    'sS: save the spectrum again',$
                    'w: click w on the bottom and top of the window to resize',$
                    'W: back to the original size']
            "RESIZE": BEGIN   ;; resize plotting window
                print,'info.action ',info.action
                IF info.action EQ "NONE" THEN BEGIN
                    info.action = "RESIZE"
                    info.tmp = cursor[1]
                ENDIF ELSE IF info.action EQ "RESIZE" THEN BEGIN
                    info.action = "NONE"
                    t1 = info.tmp
                    t2 = cursor[1]
                    IF t1 GT t2 THEN BEGIN
                        tmp=t2
                        t2=t1
                        t1=tmp
                    ENDIF
                    points.ylim=[t1,t2]
                    PlotPoints, /ylim
                ENDIF

            END
            "ORIGINALSIZE": PlotPoints   ;; resize plotting window
            "DELETE": BEGIN   ;; mask a point
                xrange = 2.2*7.172/27.*points.wpix/2.*info.pixsize
                yrange = 1.1*(points.ylim[1]-points.ylim[0])
                ;print,'range ',xrange,yrange
                dis = ((points.x-cursor[0])/xrange)^2+$
                      ((points.y-cursor[1])/yrange)^2
                mindis = min(dis,imin,/nan)
                ;print,'dis imin ',mindis,imin
                points.w[imin]=0
                points.m[imin]=0
                frame[points.fn[imin]].badpixmask[points.xn[imin],points.j+points.y1]=0
                PlotPoints, /ylim
            END
            "UNDELETE": BEGIN ;; unmask a point
                xrange = 2.2*7.172/27.*points.wpix/2.*info.pixsize
                yrange = 1.1*(points.ylim[1]-points.ylim[0])
                dis = ((points.x-cursor[0])/xrange)^2+$
                      ((points.y-cursor[1])/yrange)^2
                mindis = min(dis,imin)
                wmean = mean(points.w(where(points.w GT 0)))
                points.w[imin]=wmean
                points.m[imin]=1
                frame[points.fn[imin]].badpixmask[points.xn[imin],points.j+points.y1]=1
                PlotPoints, /ylim
            END
            "FIT": FitPoints  ;; do a new fit
            "SAVE": SaveSpectra ;; wrote spectra in ASCII files
            ELSE:
        ENDCASE
    ENDIF
ENDIF

end


PRO IrsLow, list, sl=sl, out=out, force=force, aper=aper
COMMON database
;+
; NAME:
;      IrsLow
;
; AUTHOR:
;   Dario Fadda, NHSC/Caltech, MS 100-22, Pasadena, CA 91125
;   fadda@ipac.caltech.edu
;   UPDATED VERSIONs can be found on my WEB PAGE: 
;   http://web.ipac.caltech.edu/staff/fadda/IrsLow
;
; PURPOSE: Extract spectra from low-resolution IRS Spitzer observations
;       
; CALLING SEQUENCE:
;      IrsLow, listfiles, sl=0
;
; INPUTS:
;      LIST   = Name of the list with the names of original BCDs
;      SL     = 0 for long-low, 1 for short-low spectra
;
; EXAMPLE: short-low spectra
;          IDL> list=file_search('/path/r00000000/ch0/bcd/*bcd.fits')
;          IDL> IrsLow,list,/sl
;          long-low spectra
;          IDL> list=file_search('/path/r00000000/ch2/bcd/*bcd.fits')
;          IDL> IrsLow,list
;          if files in more than one AOR
;          IDL> list=file_search('/path/r00000001/ch2/bcd/*bcd.fits')
;          IDL> list=[list,file_search('/path/r00000002/ch2/bcd/*bcd.fits')]
;          IDL> IrsLow,list
;
; PROCEDURES USED:
;       Procedures: bspline fitting from idlutils-v5_2_0
;                   in http://spectro.princeton.edu/tarballs/
;       Calibration: from spice
;       (ssc.spitzer.caltech.edu/postbcd/download-spice.html)
;
; REFERENCES: If you used the code for your analysis, please refer to
;             the paper: Fadda,D. et al., 2010, ApJ in press (arXiv:1006.2873)
;
; MODIFICATION HISTORY:
;       Written,                    Dario Fadda    v1.0    April, 2005
;       Major rewriting             Dario Fadda    v1.1     July, 2006
;       Extraction added            Dario Fadda    v1.2     July, 2006   
;       Several bugs corrected      Dario Fadda    v1.3     Sept, 2006
;       Complete rewriting          Dario Fadda    v2.0     July, 2008
;       to include new mapping mode
;       Bugs corrected, errors
;       excluding flux line (to be
;       compatible with spice)      Dario Fadda    v2.1     Oct , 2008
;       Check of flux with calibration
;       stars                       Dario Fadda    v2.2     July, 2009
;       Public release                                      June, 2010
;-
; Copyright (C) 2010, Dario Fadda
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-

if N_params() LT 1 then BEGIN
    print,'Syntax - IrsLow, listfiles [,/sl]'
    return
endif



;;; Create the info structure
channel=0
IF keyword_set(sl) THEN channel=1 
CreateInfo, list, channel

IF keyword_set(out) THEN info.out = 1
IF keyword_set(aper) THEN info.aper = 1
IF keyword_set(force) THEN info.force = 1

;;; Pixsize derived by taking the two minima in the distortion
;;; file (1st and 2nd orders) and then imposing a difference of 192"


;;; Get calibration and data
ReadData, list
ReadCalib

;;; Compute and subtract background
SubtractBackground

;;; Compute images to display for the two orders
ComputeOrder, 1
ComputeOrder, 2

;;; Launch the widgets and display orders
LaunchWidgets
DisplayImage

END
