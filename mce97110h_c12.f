	program mce97110h
! Cut apply: theta>0
! Cut at x, y_sieve to choose 1C,D,E and 2C,D,E
!________________________________________________________________________
! Check sieve IN
! Check collimator
! Check foil
! 11/14/2016
! Nguyen modified formula for scattering angle, without approximation

! Monte-Carlo of HRS hadron/right arm  spectrometer using uniform illumination.
! Nilanga Liyanage: took Alex's MC program and modified it for the broken setpum case
! please see the README file
! Since there is no magnetic model of the broken septum, the forward matrix elements
! used here from a fit to data.
! No transport through the spectrometer done, no multiple scattering effects taken into 
! account. 
!This  program has the structure the ancient SLAC monte-carlo program, uniform ( David Potterveld, Sept. 1985). 
!major changes:
!  - For UNIX and implicit none
!  - Not using the same coordinate systemes (now we are right handed)
!  - HRS setting and optics
!  - polarized 3He target setting
!  - get ride of the TopDraw setup (now PAW is used to plot output)
!  - new random numbers generator
!  - Raster option 
!  - Physics: C12 elastic cross sections
!             Radiative corrections and landau tail for elastic peak
!             3He quasi-elastic cross sections, asymmetries and external RC 
!A. Deur June 1999
!This version includes the septum magnet. March 2003
!Includes target collimators. April 2003
!Includes C12 Xsection
!________________________________________________________________________

	implicit none

	real*4 d_r,pi,r_d,res_x,res_y,res_th,res_ph,root,flagsv
C Type declarations
	character*21 qelname
	character*13 qelfile
!	logical guard

	integer i,n_trials
	integer trial,foil
!	integer stop_fid
	real rndm,rand      !random number generator must link to CERN libs
	integer*4 today(3), now(3)
	integer*4 nseed
!	real*4 delta_p,delta_t,delta_phi   	!reconstructed quantities
!	real*4 delta_x
	real*4 Etest,Etest2,eqel
	real*4 Ei,Ep,Eps,th_spec !incident energy,scat energy,spectro setting
	real*4 inl,outl ! energy losses (incoming and outgoing) in rad length 
	real*4 radw			!temp. radiation length
	real*4 rc,dE1,dE2,dE3,dE4,detot,w ! if rc=0, no elastic rad cor,
        real*4 aav,avz,zav,ang,angstrag,ioni1,ioni2  ! if rc=1, rad cor turned on
!	real*4 rcw1,rcw2,rcw3,rcw4 !weight for RC.
	real*4 cos_HRS,sin_HRS
	real*4 x,xtest,y,ytest,z,t1,t2,t3,t4,tt,u1,u2,etemp!temporaries
	real*4 t5,t6,t7,t8,et !more temporaries
	real*4 xtr,ytr,ztr			!HRS transport coords. (cm)
	real*4 xfoc,yfoc,phfoc,thfoc	!HRS coords. (cm) in focal plane
	real*4 deltat,thetat,phit,yt,zreact !HRS coords. (cm) at target (reconstructed)
	real*8 dpor,thor,phor,yor,Epor   !HRS coords. (cm) at target (original)
	real*8 wor
	real*4 dpnstrg ! delta P before straggling      
	real*4 dpp_ac				!HRS deltap/p (percent)
	real*4 dth_ac				!HRS delta theta (mr)
	real*4 dph_ac                           !HRS delta phi   (mr)
	real*4 dpp_max,dpp_min			!cut deltap/p (percent)
	real*4 dth_max,dth_min			!cut delta theta (mr)
	real*4 dph_max,dph_min			!cut delta phi   (mr)
	real*4 tgt_max,tgt_min
	real*4 w_max,w_min
	real*4 dpp				!rnd deltap/p (percent)
	real*4 dth				!rnd delta theta (mr)
	real*4 dph				!rnd delta phi   (mr)
!	real*4 z_shift				!offset of beam (cm)
!	real*4 xplane
	real*4 bcol,ccol,ecol,fcol,tanphi !variables for target collimators
	real*4 offbc,offcc,offec,offfc,ycoll !offsets collimators
	integer sieve
	real*4 dr,lc 
	real vjjl(5) ! vector for J.J LeRose forward functions.
	real x_sr_q1ex,y_sr_q1ex,x_sr_fp,y_sr_fp,p_sr_fp,t_sr_fp
	real x_sr_ep3,y_sr_ep3,x_sr_ep4,y_sr_ep4,x_sr_ep5,y_sr_ep5
	real x_sr_ep6,y_sr_ep6,x_sr_ep7,y_sr_ep7,x_sr_dent,y_sr_dent
	real x_sr_dext,y_sr_dext,x_sr_q3en,y_sr_q3en,x_sr_q3ex,y_sr_q3ex
	real yjjl(5)  ! vector for J.J LeRose reverse functions.
	real txfit,delta,theta,phi,y00 ! jjl reverse transfer functions.
	real spot_x,spot_y,tgt_l,xoff,yoff,zoff,xbeam,ybeam
	real asym_calc3,cross_section,aspin,choice,xs,asy,avasy ! variables for Xsection and asy computation
	real xom,sigma,qsq,atp,atlp,atpp,atpn,atlpp,atlpn,al,app ! variables for quasi-elastic Xsection and asy computation
	REAL XSPREV,asyprev,xomprev !temporaries for quasi-elastic interpolation
       
	real tacc,xsc,cacc,full_sa,xscount,xs_sum,cutc	! acceptance and cross section counters
	real r_phi,r_theta !angle in radian to calculate real scattering angle
	real QQ,mott
	real c_acc,c_tacc
! Function for rad corr
	real incl,intrad,ioni!,deltazero,deltazero1,deltazero2
	real xdi,xdo !average density*thickness for incoming and outgoing e-
	common/landau/ld(2,68)
	real ld ! contain the landau tail distribution
c
c NL arrays for the corners of the solid angle
c
	real y_corner(6)
	integer ii
	real ph_corner(6)
	real edge(6)
	real slope, intercept
	real t_phi_corner(5),t_theta_corner(5),t_edge(5)
	real t_xx,t_slp
	real x_sieve,y_sieve

	real tacc2,xsc2,cacc2,xs2
	real tacc3,xsc3,cacc3,xs3
	real tacc4,xsc4,cacc4,xs4
	real tacc5,xsc5,cacc5,xs5
	real tacc6,xsc6,cacc6,xs6
	real ph_min1,ph_min2,ph_min3,ph_min4,ph_min5,ph_min6
	real ph_max1,ph_max2,ph_max3,ph_max4,ph_max5,ph_max6
c
C Coefficient for forward matrix
	real coef_x(10),coef_y(12),coef_phi(10),coef_th(7)
	integer tmp_x,tmp_y,tmp_phi,tmp_th

	real dp_cut_min,dp_cut_max,yfp
	integer i_count
c	
!!!declare stuff to be able to make Raw Wise Ntuple!!

	common/pawc/blanc(150000)
	real blanc 
        real ztup(25)
        call initpaw

!zero acceptance and cross section counters
        tacc=0.
	xsc=0.
	cacc=0.
	avasy=0.

        tacc2=0.
        tacc3=0.
        tacc4=0.
        tacc5=0.
        tacc6=0.
	xsc2=0.
	xsc3=0.
	xsc4=0.
	xsc5=0.
	xsc6=0.
	cacc2=0.
	cacc3=0.
	cacc4=0.
	cacc5=0.
	cacc6=0.
C Math constants
	pi = 3.141592654
	d_r = pi/180.
	r_d = 180./pi
	root = 0.707106781 !square root of 1/2
C HRS Wire chamber resolutions (half width at one sigma.)

	res_x = 0.013	!cm
	res_y = 0.013	!cm
	res_th = 0.3	!mr
	res_ph = 0.3	!mr

C Beam position offsets 2.134 GeV, 6 degrees.
	
	xoff = -0.05043! cm
	yoff = 0.05985!cm
	zoff = -0.057 ! cm foil position offset

	open (20,file='c12.inp',status='old')
	open (25,file='cutc12.inp',status='old')
! Random seed
        call idate(today) ! select a seed for rndm
        call itime(now)
        nseed =today(1)+today(2)+now(1)+now(2)+now(3) ! 
        u1=rand(nseed)
        print*,'Print u1= ',nseed,' ',u1

! Read data lines

	read (20,*) n_trials,Ei,Ep,th_spec,dpp_ac,dth_ac,dph_ac,
     +  spot_x,spot_y,tgt_l,aspin,choice,inl,outl,xdi,xdo,rc
	read (25,*)dpp_max,dpp_min,dth_max,dth_min,dph_max,dph_min,
     +  tgt_max,tgt_min,w_max,w_min
        if (tgt_max.lt.tgt_min) then ! security
           print*, 'tgt_max smaller than tgt_min'
        endif
         if (dpp_max.lt.dpp_min) then
           print*, 'dpp_max smaller than dpp_min'
        endif
          if (dph_max.lt.dph_min) then
           print*, 'dph_max smaller than dph_min'
        endif

 	if (dth_max.lt.dth_min) then
           print*, 'dth_max smaller than dth_min'
        endif
	if (tgt_max.gt.tgt_l/2.) then
	   print*, 'target is not completely illuminated'
	endif
	if (abs(tgt_min).gt.tgt_l/2.) then
	   print*, 'target is not completely illuminated'
	endif
	if (dpp_max.gt.dpp_ac/2.) then
	   print*, 'momentum range is not completely illuminated'
	endif
	if (abs(dpp_min).gt.dpp_ac/2.) then
	   print*, 'momentum range is not completely illuminated'
	endif
	if (dph_max.gt.dph_ac/2.) then
	   print*, 'phi range is not completely illuminated'
	endif
	if (abs(dph_min).gt.dph_ac/2.) then
	   print*, 'phi range is not completely illuminated'
	endif
	if (dth_max.gt.dth_ac/2.) then
	   print*, 'theta range is not completely illuminated'
	endif
	if (abs(dth_min).gt.dth_ac/2.) then
	   print*, 'theta range is not completely illuminated'
	endif			! end of security checks

! Read coefficient for forward matrix
	open(30,file='infiles/x_coef.inp',status='old')
	open(31,file='infiles/y_coef.inp',status='old')
	open(32,file='infiles/phi_coef.inp',status='old')
	open(33,file='infiles/th_coef.inp',status='old')
	do tmp_x = 1,10
	   read(30,*) coef_x(tmp_x)
	enddo
	do tmp_y = 1,12
	   read(31,*) coef_y(tmp_y)
	enddo
	do tmp_phi = 1,10
	   read(32,*) coef_phi(tmp_phi)
	enddo
	do tmp_th = 1,7
	   read(33,*) coef_th(tmp_th)
	enddo

! End read coefficent	   
 	
	etemp=Ei

	
	th_spec = th_spec*d_r ! in rad
	cos_HRS = cos(th_spec)
	sin_HRS = sin(th_spec)
	spot_x=spot_x*100. ! in cm.
	spot_y=spot_y*100. ! in cm.
	
	!Keep cuts on target length (z) and not on y_target 12/20/2005
	!tgt_min=tgt_min*abs(sin_HRS)!tgt_min is now the cut on y_target (HRS
	!tgt_max=tgt_max*abs(sin_HRS)!frame) not a cut on target length (z) anymore 
	if (rc.eq.1.) then ! fill the common block for landau distribution
	   call stuffit
	endif
C look at physics database
****elastic is in a subroutine *************************
***quasi-elastic****************************************
**this is now not valid anymore***************
	if (choice.eq.2) then	!quasi-elastic
	   etest2=0.
	   open(23,file='qel.dir/qel_en.dat',status='old')
	do i=1,100 ! chose the right incident energy file
 	   read(23,*,end=998) Etest,qelfile
	   if (abs(etest-Ei).le.abs(etest2-Ei)) then
	   eqel=Etest
	   qelname= 'qel.dir/' //  qelfile
	   endif
	   etest2=etest
 998	continue
	enddo

	open(24,file=qelname,status='old')!contains the quasi-elastic data
	!print*,qelname 
	endif
***end of quasi-elastic****************************************
	if (choice.eq.0) then   !phase-space
c
c Calculate the full solid angle for the box we have defined
c dp values are given as percent, need to be divided by 100
c theta and phi are given as mrad, need to be divided by 1000.	   
c
	   full_sa = (Ep*dpp_max-Ep*dpp_min)/100.0
     +               *(dth_max-dth_min)/1000.0
     +               *(dph_max-dph_min)/1000.0     
c
	   write(6,*) "Full solid angle is ",full_sa," Gev. st_rad"
	endif
c
	if (choice.eq.1) then   !phase-space for elastic
c
c Calculate the full solid angle for the box we have defined
c	   
	   full_sa = (dth_max-dth_min)/1000.0
     +               *(dph_max-dph_min)/1000.0     
c
	   write(6,*) "Full solid angle is ",full_sa," st_rad"
	endif

c
C------------------------------------------------------------------------C
C                       Top of Monte-Carlo loop                          C
C------------------------------------------------------------------------C

	do foil = 0,0
	do trial = 1,n_trials

C Pick two indep., normal dist. numbers, used only if raster is off.
 56	u1=rand()!rndm(trial)
	u2=rand()!rndm(trial+1)
!	t1=sqrt(-2.*log(u1))*cos(2*pi*u2) !t1 norm dist (cf notes Fonvieille on  MC)
!	u1=rndm(trial+n_trials) 
!	u2=rndm(n_trials+trial+1)
!	t2=sqrt(-2*log(u1))*cos(2*pi*u2)

C Pick starting point within target. Z is picked uniformly, X and Y are
C chosen as truncated gaussians. (Truncated to 3 sigma) if raster off
C If CIRCULAR raster is on, uniform dist. 
C Units are cm.
     
!	x = t1*spot_x+xoff	! if raster off, unit in cm
!	y = t2*spot_y+yoff	! if raster off, in cm

	if (((2*u1-1)**2.+(2*u2-1)**2.).gt.1.) then! if circular raster on (uniform dist), 
                   ! in cm. for a squared uniform dist., comment the 3 lines 
	goto 56
	endif

	x=(u1-0.5)*spot_x+xoff !*2.
 	y=(u2-0.5)*spot_y+yoff !*2.
	xbeam=x
	ybeam=y

	z=rand()!rndm(trial+2)
	z=(z-0.5)*tgt_l+(zoff+(foil-1)*10.0) ! in cm 
	
C       Transform from target to HRS transport coordinate 

!	   xtr = -y
!	   ytr = x*cos_HRS - z * sin_HRS ! changed -x to +x 12/07/2005
	   !ztr = z * cos_HRS - x * sin_HRS

C Pick scattering angles and DP/P after to have computed energy losses
C picks a gaussian dp/p center around dp/p_scat
C pick a gaussian for Beam energy dispertion
	    Ei=Etemp !last Ei was shift by energy loss and dispertion. Come back at the set one 
C********************
C This is where problem with radiation correction come from
C If I use rand(), it will create negative energy
C If I use rndm, no problem 
	    u1=rndm(trial+10)
	    u2=rndm(trial+11)
  	    t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist

	    Ei=t1*3*10.**(-5.)*Ei+Ei ! beam dispertion is taken as 3.e-5
	    Et=Ei		! used in W - M calculation
	    
	    
	    dth =rand()!rndm(n_trials+4)
	    dth =(dth-0.5)*dth_ac

	    dph =rand()!rndm(n_trials+5) 
            dph =(dph-0.5)*dph_ac



	    xtr = -y -dth/1000.0*(z*cos_HRS/cos(atan(dph/1000.0)))!transport coords
	    ytr=sin(th_spec+atan(dph/1000.0))/cos(atan(dph/1000.0))*
     >           (x/tan(th_spec+atan(dph/1000.0))-z)
!	    ang=th_spec+dph/1000.-(dth/1000.)**2.! now second order theta
! Without approximation
	    r_phi=dph/1000.0
	    r_theta =dth/1000.0
	ang=acos(cos(r_theta)*(cos_HRS*cos(r_phi)+abs(sin_HRS)
     +    *sin(r_phi))
     +  /sqrt(sin(r_phi)*sin(r_phi)*cos(r_theta)*cos(r_theta)
     +  *cos_HRS*cos_HRS+sin(r_theta)*sin(r_theta)*cos(r_phi)*cos(r_phi)
     + *cos_HRS*cos_HRS+cos(r_phi)*cos(r_phi)*cos(r_theta)*cos(r_theta))
     +	)

! With approximation in small angle
!	    ang=acos((cos_HRS + abs(sin_HRS)*abs(dph/1000.))
!     +          /sqrt(1+(dph/1000.)**2.+(dth/1000.)**2.))
	      if (rc.eq.1.) then ! Rad cor enable. pick up an
		 !!energy loss on incoming particle!!
	         ! pick up an dEi energy loss du to external 
		 aav=4. ! average A for incoming e-
		 zav=2.  ! average Z for incoming e-
		 dE1=incl(Ei,zav,inl)  
		 Ei=Ei-dE1 
	         ! pick up an dEi energy loss du to internal 	 	 
 		 dE2=intrad(Ei,abs(ang))
		 ioni1=ioni(Ei,aav,zav,xdi)! ioni is loss by ionization 
                 !deltazero1=deltazero(Ei,aav,zav,xdi)!most probab loss by ioni
		 Ei=Ei-dE2-ioni1 !-0.000207

                 !now chose or compute scattered particle energy
		 if (choice.eq.1) then ! elastic
		 !have to compute Eps with new Ei
	         Eps=Ei/(1+2*Ei/11.178*(sin(ang/2.))**2.)!11.178 GeV mass c12
	         Epor=Etemp/(1+2*Etemp/11.178*(sin(ang/2.))**2.)!without RC and resolution
		  else ! quasi elastic or phase space 
	          dpp =rand()!rndm(n_trials+3)
	          Eps =Ep*(1+(dpp-0.5)*dpp_ac/100.)
	          Epor =Eps !without RC and resolution
		 endif   
		 !!!!!!!

		 ! now outgoing e- energy loss due to internal
 		 dE3=intrad(Eps,abs(ang)) 
		 Eps=Eps-dE3	
		! external loss of outgoing e-
		 aav=9.14 ! average A for outgoing e-
		 zav=4.77  ! average Z for outgoing e-
	    if(Eps<0) print*,'*************= ',Eps	

		 dE4=incl(Eps,zav,outl)
		 ioni2=ioni(Eps,aav,zav,xdo)
                 !deltazero2=deltazero(Eps,aav,zav,xdo)
		 Eps=Eps-dE4-ioni2 ! -0.00007
	         ! Now compute the final relative momentum dpp with resolution
!                  Will correct for HRS mom resolution when vdc resolution is 
!                  taken care of on lines 650-676
	         u1=rand()!rndm(trial+16)
	         u2=rand()!rndm(trial+17)
  	         t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist
	         t1=t1*6.*10.**(-4.)*Ep !relative HRS mom res taken as 3e-4
		 Eps=Eps+t1
		 dpp=(Eps-Ep)*100./Ep	!   dpp has to be in % 
	                   !   (converted in fractional value in 
	                   !    subroutine scavec for JJL function)
		  !print*, de1,de2,de3,dpp

	!write energy loss and equialent w to file. Effect of the incoming 
        !energy loss is corrected for recoil.
		 dEtot=(de1+de2)/(1+2*Ei/11.17808*(sin(ang/2.))**2.)+de3+de4
     ++ioni1+ioni2

	    else ! rad Cor not on
	       u1=rand()!rndm(trial+16)
	       u2=rand()!rndm(trial+17)
	       t1=sqrt(-2.*log(u1))*cos(2.*pi*u2) ! gaussian dist
	       t1=t1*7.*10.**(-4.)*Ep !relative HRS mom res taken as 3e-4
	       if (choice.eq.1) then ! elastic
	        Epor=Etemp/(1+2*Etemp/11.178*(sin(ang/2.))**2.)!without RC and resolution 
	        Eps=Ei/(1+2*Ei/11.178*(sin(ang/2.))**2.)!11.178 GeV mass c12
	        Eps=Eps+t1
	        dpp=(Eps-Ep)*100./Ep 
	        else ! quasi elastic or phase space 
	        dpp =rand()!rndm(n_trials+3)
	        Eps =Ep*(1+(dpp-0.5)*dpp_ac/100.)
	        Epor =Eps!without RC and resolution
	        dpp=(Eps-Ep)*100./Ep 
	       endif 

	    endif ! end of energy loss computation due to RC


C origin coordinates in the HRS frame

	   yor=ytr
	   phor=dph
	   thor=dth
	   dpor=(Epor-Ep)*100./Ep ! otherwise, shifted by radiative proces

	   dpnstrg=dpp
	   wor=sqrt(11.178**2+2*11.178*(Etemp-Epor)
     +	   -4*Etemp*Epor*(sin(-1*ang/2.))**2.)-11.178
C
C  Calculate reactz from reconstructed quantities at target
	  zreact = -(yor/100.0)*cos(atan(phor/1000.0))/(sin(th_spec + 
     +    atan(phor/1000.0)))+xbeam/(100.0*tan(th_spec 
     +    + atan(phor/1000.0)))


C Project trajectory to sieve slit
C sieve slit box to sieve slit location

	tt = 80.06		! cm 
	call project(xtr,ytr,tt,dth,dph)

	   sieve=1
	   if (sieve.eq.1) then ! sieve slit is in
	      call sieve_slit(xtr/100.0,ytr/100.0,flagsv)
	      if (flagsv.lt.1.0) then
		 goto 500
	      endif
	   endif

!sieve thickness
	tt = 0.5		! cm 
	call project(xtr,ytr,tt,dth,dph)

C	   sieve=1
C	   if (sieve.eq.1) then ! sieve slit is in
	      call sieve_slit(xtr/100.0,ytr/100.0,flagsv)
	      if (flagsv.lt.1.0) then
		 goto 500
	      endif
C	   endif

C Cut on theta target
	   if(thor.le.10.or.thor.ge.70) then
	      goto 500
	   endif

C Cut on momentum range for phase space: -3.5 to 4.5
        if(dpor.le.-3.5.or.dpor.ge.4.5) then
	   goto 500
	endif
C Cut on z_target: +/- 18cm
        if(z.le.-18.0.or.z.ge.18.0) then
	   goto 500
	endif
c use this cut for foil target runs only, 3 12C foils at +/- 10 cm and 0
c
c	  if ((abs(zreact-0.1).ge.0.01).and.(abs(zreact+0.1).ge.0.01)
c     +         .and.(abs(zreact).ge.0.01))goto 500

c NL: Apply box cuts given in the cuts file for zreact,th and ph

c Ton: make sure use correct collimator for each case
* optic: sieve collimator: only cut for optic
* acceptance: target collimator
c Apply 6 degrees collimator block cuts
c
	   if ((abs(th_spec).ge.(5.75*d_r)).and.
     @     (abs(th_spec).le.(6.25*d_r))) then    
!!!!!!!!ice cone collimator!!!!!!!!!
!C check if it going through the middle collimator
	   ecol=18.99 !position collimator along the beam line
	   fcol=2.326 !dist from collimator to beam axis
	    tanphi=(fcol+x)/abs(z-ecol)/cos(thor/1000.)
C	     if((tan(abs(th_spec+phor/1000.))).ge.tanphi) goto 500

!C check if it going through the downstream window collimator
	   bcol=33.96 !z position of collimator
	   ccol=2.85 !dist top of collimator to beam axis
	    tanphi=(ccol+x)/abs(bcol-z)/cos(thor/1000.)
C	    if((tan(abs(th_spec+phor/1000.))).le.tanphi) goto 500
c	    endif !end of 6 degree 

!!!!!!!!standard target!!!!!!!!!!!!
!C check if it going through the first edge of middle collimator
c ton	    ecol=0.0293		!position collimator along the beam line
c ton	    fcol=-0.0107	!dist from collimator to beam axis
c ton	    ycoll=yor/100.+ecol*phor/1000.
c ton	    if (ycoll.le.fcol) goto 500
C       check if it going through the second edge of middle collimator
c ton	    ecol=0.0594		!position collimator along the beam line
c ton	    fcol=-0.0117	!dist from collimator to beam axis
c ton	    ycoll=yor/100.+ecol*phor/1000.
c ton	    if (ycoll.le.fcol) goto 500
c
!C check if it going through the first edge of downstream window collimator
c ton	    bcol=0.4197		!z position of collimator
c ton	    ccol=0.00692	!dist top of collimator to beam axis
c ton	    ycoll=yor/100.+bcol*phor/1000.
c ton	    if (ycoll.ge.ccol) goto 500
C check if it going through second edge of the downstream window collimator
c ton	    bcol=0.450		!z position of collimator
c ton	    ccol=0.0076  		!dist top of collimator to beam axis
c ton	    ycoll=yor/100.+bcol*phor/1000.
c ton	    if (ycoll.ge.ccol) goto 500
c
c NL: for the optics runs use a wide cut at the location of the collimator
c 0.025 at 0.5 m  corresponds to roughly +/- 10 cm target and 30 msr phi
c
	      bcol = 0.45
	      ccol = 0.05
	      ycoll=yor/100.+bcol*phor/1000.
c
	      if (abs(ycoll).ge.ccol)then
		 goto 500
	      endif
	 endif			!end of 6 degree 



C
C Now calculate the real focal plane quantities for the broken septum
C case using NL forward matrix from target to focal plane
C
C	  convert the target variables  to correct units first

	  phit=phor/1000.0
	  thetat=thor/1000.0
	  yt = yor/100.0
	  deltat=dpp/100.0
! 97 data and 39 parameters	  
c============== y focal ========================
          yfoc =  coef_y(1)
     +           +coef_y(2)*deltat
     +           +coef_y(3)*phit
     +           +coef_y(4)*phit*deltat
     +           +coef_y(5)*phit*phit
     +           +coef_y(6)*yt
     +           +coef_y(7)*yt*deltat
     +           +coef_y(8)*phit*yt
     +           +coef_y(9)*thetat
     +           +coef_y(10)*thetat*thetat
     +           +coef_y(11)*thetat*deltat
     +           +coef_y(12)*thetat*phit
   
c============== x focal ========================
          xfoc =  coef_x(1)
     +           +coef_x(2)*deltat
     +           +coef_x(3)*deltat*deltat
     +           +coef_x(4)*phit
     +           +coef_x(5)*phit*phit
     +           +coef_x(6)*yt
     +           +coef_x(7)*thetat
     +           +coef_x(8)*thetat*thetat
     +           +coef_x(9)*thetat*yt
     +           +coef_x(10)*thetat*phit
c============== theta focal ========================
          thfoc = coef_th(1)
     +           +coef_th(2)*deltat
     +           +coef_th(3)*phit
C     +           +coef_th(4)*phit*phit
     +           +coef_th(4)*yt
     +           +coef_th(5)*thetat
     +           +coef_th(6)*thetat*yt
     +           +coef_th(7)*thetat*phit                        
c==============phi focal ==========
          phfoc = coef_phi(1)
     +           +coef_phi(2)*deltat
     +           +coef_phi(3)*phit
     +           +coef_phi(4)*phit*deltat
     +           +coef_phi(5)*phit*phit
     +           +coef_phi(6)*yt
     +           +coef_phi(7)*yt*deltat
     +           +coef_phi(8)*thetat
     +           +coef_phi(9)*thetat*deltat
     +           +coef_phi(10)*thetat*phit

c	  
C Reconstructed ytg, thetatg, phitg

C	yt = 0.00075+1.36*yfoc-1.59*phfoc+0.0116*xfoc+2.28*thfoc
C     + +6.77*yfoc*phfoc-39.22*yfoc*thfoc-67.98*phfoc*thfoc
C     + +6.22*thfoc*thfoc+28.01*phfoc*phfoc-13.36*yfoc*yfoc

C	thetat =0.037990-0.004055*xfoc+0.849883*yfoc+0.07456*thfoc
C     + -0.163464*phfoc
 
C	phit=0.001692-0.012293*xfoc-0.007802*yfoc-0.40254*thfoc
C     + +0.627065*phfoc-5.927*yfoc*phfoc-6.172*yfoc*thfoc
C     + -11.326*yfoc*yfoc+6.111*phfoc*phfoc

C	deltat=0.004664+0.075442*xfoc+0.098884*yfoc+0.51163*thfoc
C     + +0.389*thfoc*xfoc+5.7609*thfoc*thfoc

! NT: new corner (based on Nilanga method)
c$$$	y_corner(1) = -0.187788*deltat + 0.007736
c$$$	y_corner(2) = -0.182372*deltat - 0.013877
c$$$	y_corner(3) =  0.012436*deltat - 0.015549
c$$$	y_corner(4) =  0.514038*deltat + 0.013196
c$$$	y_corner(5) = -0.259872*deltat + 0.026024
c$$$	y_corner(6) = -0.453814*deltat + 0.019524
c$$$	
c$$$	ph_corner(1) = -0.088333*deltat + 0.023010
c$$$	ph_corner(2) = -0.479423*deltat + 0.038124
c$$$	ph_corner(3) =  0.201923*deltat - 0.016102
c$$$	ph_corner(4) =  0.629295*deltat - 0.030005
c$$$	ph_corner(5) = -0.015769*deltat + 0.045560
c$$$	ph_corner(6) = -0.064167*deltat + 0.047718

	y_corner(1) = -0.227484*deltat + 0.00965004
	y_corner(2) = -0.162115*deltat - 0.017178
	y_corner(3) =  0.203004*deltat - 0.0157547
	y_corner(4) =  0.524234*deltat + 0.0129056
	y_corner(5) = -0.260981*deltat + 0.024819
	y_corner(6) = -0.394433*deltat + 0.0205985
	
	ph_corner(1) = -0.120625*deltat + 0.0271181
	ph_corner(2) = -0.0885365*deltat + 0.012184
	ph_corner(3) = 0.209205*deltat - 0.0203365
	ph_corner(4) = 0.677652*deltat - 0.0290798
	ph_corner(5) = 0.105192*deltat + 0.0445469
	ph_corner(6) = -0.0579263*deltat + 0.0470789

          do i=1,6
             if((i.eq.2).or.(i.eq.4)) then
                slope = (y_corner(i+1)-y_corner(i))/
     +                  (ph_corner(i+1)-ph_corner(i))
                intercept = y_corner(i) - slope*ph_corner(i)
                edge(i) = slope*phfoc + intercept
             else
                if(i.eq.6) then
                   slope = (ph_corner(1)-ph_corner(i))/
     +                  (y_corner(1)-y_corner(i))
C                   intercept = ph_corner(i) - slope*y_corner(i)
                else
                   slope = (ph_corner(i+1)-ph_corner(i))/
     +                  (y_corner(i+1)-y_corner(i))
C                   intercept = ph_corner(i) - slope*y_corner(i)
                endif
		intercept = ph_corner(i) - slope*y_corner(i)
                edge(i) = slope*yfoc + intercept
             endif
          enddo


c
c Now apply the solid angle cut


	if((phfoc.gt.edge(1).and.phfoc.gt.edge(6)).or.
     + phfoc.lt.edge(3).or.phfoc.gt.edge(5).or.
     + yfoc.lt.edge(2).or.yfoc.gt.edge(4)) then
	   goto 500
	endif



C End reconstructed
c convert the units back to the ones used in this program 
c
	  phit=phit*1000.0
	  thetat=thetat*1000.0
	  yt = yt*100.0
	  deltat=deltat*100.0
c
C Monte-Carlo trial is a 'GOOD' event

! Compute theoritical cross section and assymetries for each event
! kinematic Need to use dp/p before outgoing energy loss for Xsec computation
!	angstrag=th_spec+phor/1000.+(thor/1000.)**2. 

! Angle formula without approximation
	angstrag=acos(cos(thetat/1000.)*(cos_HRS*cos(phit/1000.)
     +  +abs(sin_HRS)*sin(phit/1000.))/sqrt(sin(phit/1000.)
     +  *sin(phit/1000.)*cos(thetat/1000.)*cos(thetat/1000.)
     + *cos_HRS*cos_HRS+sin(thetat/1000.)*sin(thetat/1000.)
     + *cos(phit/1000.)*cos(phit/1000.)*cos_HRS*cos_HRS
     + +cos(phit/1000.)*cos(phit/1000.)*cos(thetat/1000.)
     + *cos(thetat/1000.)))

! With approximation in calculate angle
!	  angstrag=acos((cos_HRS + abs(sin_HRS)*phit/1000.)
!     +           /sqrt(1+(phit/1000.)**2.+(thetat/1000.)**2.))
	  Eps=Ep*(1+dpnstrg/100.)
	  w=sqrt(11.178**2+2*11.178*(Etemp-Ep*(1+dpp/100.))
     +    -4*Etemp*Ep*(1+dpp/100.)*(sin(-1*angstrag/2.))**2.)-11.178

          QQ=4*Ei*Eps*(sin(-1*ang/2.))**2
	  mott=(6.*1./137.*cos(ang/2)/2./Ei/(sin(ang/2))**2.)**2.	

!============================================
	  ph_max1 = 23
	  ph_min1 = 17

	  ph_max2 = 17
	  ph_min2 = 11

	  ph_max3 = 11
	  ph_min3 = 5

	  ph_max4 = 5
	  ph_min4 = 0

	  ph_max5 = 0
	  ph_min5 = -8

	  ph_max6 = 4
	  ph_min6 = -3

	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max1
     +  .and.phor.ge.ph_min1.and.w.le.0.01)then
	     cacc=cacc+1
	  endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max2
     +  .and.phor.ge.ph_min2.and.w.le.0.01)then
	     cacc2=cacc2+1
	  endif

	  
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max3
     +  .and.phor.ge.ph_min3.and.w.le.0.01)then
	     cacc3=cacc3+1
	  endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max4
     +  .and.phor.ge.ph_min4.and.w.le.0.01)then
	     cacc4=cacc4+1
	  endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max5
     +  .and.phor.ge.ph_min5.and.w.le.0.01)then
	     cacc5=cacc5+1
	  endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max6
     +  .and.phor.ge.ph_min6.and.w.le.0.01)then
	     cacc6=cacc6+1
	  endif
C Count tacc to calculation cross section
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max1
     +  .and.phor.ge.ph_min1.and.wor.le.0.01)then
           tacc=tacc+1
	endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max2
     +  .and.phor.ge.ph_min2.and.wor.le.0.01)then
           tacc2=tacc2+1
	endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max3
     +  .and.phor.ge.ph_min3.and.wor.le.0.01)then
           tacc3=tacc3+1
	endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max4
     +  .and.phor.ge.ph_min4.and.wor.le.0.01)then
           tacc4=tacc4+1
	endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max5
     +  .and.phor.ge.ph_min5.and.wor.le.0.01)then
           tacc5=tacc5+1
	endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max6
     +  .and.phor.ge.ph_min6.and.wor.le.0.01)then
           tacc6=tacc6+1
	endif
! increment tot acceptanace counter if the original quantities
! are in the chosen cuts
***********elastic**********************************
	if (choice.eq.1) then ! compute elastic
	   xs=cross_section(Ei,abs(ang))
	   if(rc.eq.1)then	!correct for vaccum polar and vertex 
	      xs=xs*(1+0.004647*(13./12.*log(QQ/(0.000511)**2.)-14./9.))
	   endif
	   asy=0
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max1
     +  .and.phor.ge.ph_min1.and.w.le.0.01)then
	      xsc=xs+xsc
              avasy=xs*asy+avasy
	   endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max2
     +  .and.phor.ge.ph_min2.and.w.le.0.01)then
	      xsc2=xs+xsc2
	   endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max3
     +  .and.phor.ge.ph_min3.and.w.le.0.01)then
	      xsc3=xs+xsc3
	   endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max4
     +  .and.phor.ge.ph_min4.and.w.le.0.01)then
	      xsc4=xs+xsc4
	   endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max5
     +  .and.phor.ge.ph_min5.and.w.le.0.01)then
	      xsc5=xs+xsc5
	   endif
	 if(thor.le.dth_max.and.thor.ge.dth_min.and.phor.le.ph_max6
     +  .and.phor.ge.ph_min6.and.w.le.0.01)then
	      xsc6=xs+xsc6
	   endif	     
!	print*, 'XS=',xs,' Ei=',Ei,' Theta=',thetat
!	print*, 'asy=',asy
**********quasi-elastic*****************************	
	else 
	   if (choice.eq.2) then ! compute quasi-elastic
	      if((Ei/(1.+2.*Ei/11.178*(sin(ang/2.))**2.)-eps).le.0.0054)
     +              then 
				! below quasi-elastic threshold
		                ! no need to compute it		 
		 xs=0.0
		 asy=0.
		 !print*,'trop petit'
	      else ! compute cross section and asymetry

              ! initialyze values for interpolation.
		 rewind(24)	! come back at the beginning of the file
		 read(24,*)xom,sigma,qsq,atp,atlp,atpp,
     +           atpn,atlpp,atlpn,al,app
		 xomprev=0.	!we don't want it to be xom as initial value !
		 if (aspin.eq.0.) then ! spin orientation
		    asyprev=al	! longitudinal asymetries
		    xsprev=sigma
		 else
		    asyprev=app !transverse
		    xsprev=sigma
		 endif
	      ! end of initialization

		 rewind(24)	! come back at the beginning of the file
		 do i=1,500
		    read(24,*,end=997)xom,sigma,qsq,atp,atlp,atpp,
     +              atpn,atlpp,atlpn,al,app
		    xom=xom/1000. ! now xom in GeV (it was in MeV in the file)
	                    ! xom ordered in file by increasing values
		    if(xom-(Ei-Eps).ge.0)then ! take closest sup energy lost value 
				! and go out of the loop 
		 !print*, 'going out'
		       if (aspin.eq.0.) then ! spin orientation
			  asy=al ! longitudinal asymetries
			  xs=sigma
		       else
			  if(aspin.eq.90.) then
			     asy=app !transverse
			     xs=sigma
			  else 
			     print*,'no good angle for q-elastic (must be 0 or 90)'
			     goto 996
			  endif 
		       endif
		       goto 995
		    endif
		    xomprev=xom
		    if (aspin.eq.0.) then ! spin orientation
		       asyprev=al ! longitudinal asymetries
		       xsprev=sigma
		    else
		       asyprev=app !transverse
		       xsprev=sigma
		    endif
 997		    continue
		 enddo

 995		 continue	  
	    !print*, 'qtties:',asy,xs,xom
		 xs=(xs-xsprev)/(xom-xomprev)*(Ei-Eps-xomprev)+xsprev !interpolate
		 asy=(asy-asyprev)/(xom-xomprev)*(Ei-Eps-xomprev)+asyprev
	    !print*,'prev qtties:',asyprev,xsprev,xomprev
	      endif
 	    !print*, Ei-Eps,xom,asy,xs  
	   else   
*********phase space*********************************************
	      xs=0.
	      asy=0.	   
	   endif !if elastic
	endif	!if quasi-elastic

C Monte-Carlo is finished, copy result to Raw Wise Ntuple
	ztup(1)=zreact
	ztup(2)=xs
	ztup(3)=asy
	ztup(4)=w
	ztup(5)=mott
	ztup(6)=yt
	ztup(7)=phit
	ztup(8)=thetat
	ztup(9)=deltat
	ztup(10)=xfoc
	ztup(11)=yfoc
	ztup(12)=phfoc
	ztup(13)=thfoc
	ztup(14)=Eps
	ztup(15)=Epor
	ztup(16)=Etemp
	ztup(17)=ang
	ztup(18)=x_sieve
	ztup(19)=y_sieve
	ztup(20)=angstrag
	ztup(21)=yor
	ztup(22)=wor
	ztup(23)=dpor
	ztup(24)=thor
	ztup(25)=phor
	call hfn(1,ztup)
 500	continue


c
	enddo  !end of trials loop


	enddo  !end of foils loop
	call endpaw

	close(20)
	close(21)
	close(22)
	close(23)
	close(24)
	close(25)
	close(30)
	close(31)
	close(32)
	close(33)
C output physics choice.
	if (choice.eq.0) then
	   print*, 'phase space'
	   print*,'acceptance in the chosen cut: ', cacc/tacc
	   print*,'total events in the full SA: ', tacc
	   print*,'Events that are within our SA: ', cacc
	   print*,'Available SA: ', cacc/tacc*full_sa,' Gev st.rad'
	endif
	if (choice.eq.1) then
	   print*, 'elastic'
	   print*,'Xsection= 1',xsc/tacc,' microbarn'
	   print*,'Xsection= 2',xsc2/tacc2,' microbarn'
	   print*,'Xsection= 3',xsc3/tacc3,' microbarn'
	   print*,'Xsection= 4',xsc4/tacc4,' microbarn'
	   print*,'Xsection= 5',xsc5/tacc5,' microbarn'
	   print*,'Xsection= 6',xsc6/tacc6,' microbarn'
	   print*,'acceptance in the chosen cut: ', cacc/tacc
	   print*,'total events in the full SA: ', tacc
	   print*,'Events that are within our SA: ', cacc
	   print*,'Available SA: ', cacc/tacc*full_sa,' st.rad'
	   print*,'number of counts in acceptance ',xsc
	endif
	if (rc.eq.1) then
	      print*, 'radiative  correction on'
	      else
	      print*, 'radiative  correction off'	 
	endif
	if (choice.eq.2) then
	   print*, 'quasi-elastic'
	endif
	if (sieve.eq.1) print*, 'sieve slit in'
 996	continue
      	end





*******************************************************************************

        subroutine musc(rad_len,dth,dph,Ep,th_spec,z)
C+_____________________________________________________________________
!
! MUSC - Simulate multiple scattering of an electron in the HRS
!   spectrometer. This subroutine is a repaired and simplified version
!   of the old SLAC subroutine (used in MONTP) by the same name.
!   This subroutine will only calculate multiple scattering for electrons
!   in the HRS spectrometer.
!
! ASSUMPTIONS: DTH and DPH given in milli-radians, RAD_LEN in radiation
!   lengths. The formula used is due to Rossi and Greisen (See the book
!   by Segre, NUCLEI AND PARTICLES, 1982, p. 48.) The formula assumes a
!   gaussian distribution for the scattering angle. It is further assumed
!   that the angles DTH and DPH are the delta-theta and delta-phi angles
!   of an electron in the HRS, located at an angle with respect to
!   the beam direction that is large enough so that DTH-DPH space can
!   be approximated as a cartesian space.
!   Again, not the same subroutine as in uniform since we are in the HRS right
!   handed frame, not in the uniform left handed frame!
! from D. Potterveld 
C-_____________________________________________________________________

        implicit none
	real es,rndm,t1,pi,phi_sigma,phi_scat,theta_scat,u1,u2,rand
        real*4 rad_len,dth,dph,z
        real*4 Ep,th_spec
        pi = 3.141592654
        es = 13.6E-03 ! so energy has to be GeV

C Compute scattering angles, PHI_SCAT from a gaussian distribution,
C THETA_SCAT from uniform distribution.
        u1=rand()!rndm(1)
	u2=rand()!rndm(2)
	t1=sqrt(-2*log(u1))*cos(2*pi*u2)
C       theta_sigma = es*sqrt(rad_len)/Ep
C    CHANGED TO HAVE THE CORRECTION TERM ( FROM THE BOOKLET) 2/13/89
        phi_sigma = es*sqrt(rad_len)/Ep*z*(1.+0.038*LOG(RAD_LEN))
        phi_scat = t1*phi_sigma/2.  
        theta_scat = 2*pi*rand()!rndm(3)

C Compute new trajectory angles (units are mr)
        dth = dth + phi_scat*cos(theta_scat)*1000.
        dph = dph + phi_scat*sin(theta_scat)*1000.
        return
        end
*****************************************************************************

       subroutine project(xtr,ytr,z_drift,dth,dph)

C+___________________________________________________________________________
!
!   Calculate new HRS coordinates after drifting in a field
!   free region for a distance z_drift. It is assumed that all X, and Y
!   variables have the same (linear) dimensions, and that the angles DTH and
!   DPH are in units of milliradians and are close enough to zero that a
!   paraxial approximation may be made.
!
!   not the same subroutine as in uniform since we are in the HRS right
!   handed frame, not in the uniform left handed frame
C-___________________________________________________________________________

        implicit  none

        real*4 z_drift,xtr,ytr,dth,dph

        xtr = xtr + (dth/1000.) * z_drift
        ytr = ytr + (dph/1000.) * z_drift

        return
        end
*******************************************************************************
c
	subroutine sieve_slit(x,y,flag)
c
c added 050305,  cuts on 49 sieve holes (thin sieve)
c v.sulkosky 050405
c cut rays at the sieve slit to simulate sieve slit data
c
c x, xsieve
c y, ysieve
	real*4 x
	real*4 y
	real*4 flag
	real*4 ans
	integer i,j

c     
c     r-arm sieve slit:
c Old survey A870
c	real ysvpos(7)/0.0151401,0.0103649,0.0055897,0.0008145,
c     >  -0.0053069,-0.0114283,-0.0175497/ ! in m
c New survey A870r
	real xsvpos(7)/0.0380788,0.0247692,0.0114596,-0.00185,
     >  -0.0151596,-0.0284692,-0.0417788/ ! in m
	real ysvpos(7)/0.0155925,0.0108173,0.0060421,0.0012669,
     >  -0.0048545,-0.0109759,-0.0170973/ ! in m
c without offset
c	real xsvpos(7)/0.0399288,0.0266192,0.0133096,0.0,
c     >  -0.0133096,-0.0266192,-0.0399288/ ! in m
c	real ysvpos(7)/0.0143256,0.0095504,0.0047752,0.0,
c     >  -0.0061219,-0.0122438,-0.0183657/ ! in m

	real holesize(2)/0.000698,0.001346/ ! radius of hole in m
      
	i=1
	j=1
	flag=0.0
	do while ((i.lt.8) .and. (flag.eq.0.0)) ! sieve row
	   do while ((j.lt.8) .and. (flag.eq.0.0)) ! sieve column
				! check for large holes
	      if((i.eq.4 .and. j.eq.4) .or. (i.eq.6 .and. j.eq.5)) then
		 ans=sqrt((x-xsvpos(i))**2+(y-ysvpos(j))**2)
		 if(ans.gt.holesize(2)) then
		    flag = 0.0
		 else 
		    flag = 1.0
		 endif
	      else
		 ans=sqrt((x-xsvpos(i))**2+(y-ysvpos(j))**2)
		 if(ans.gt.holesize(1)) then
		    flag = 0.0
		 else 
		    flag = 1.0
c       write(*,*)' ans = ',ans, ' flag = ',flag
		 endif
	      endif		! sieve row loop
	      j=j+1
	   enddo		! sieve column loop
	   i=i+1
	   j=1
	enddo
	
	return
	end
c
c       end of sieve slit

*******************************************************************************
!!!!!!!!function to compute C12 elastic cross section
	real function cross_section(es,theta)

	implicit none
C     ELASTIC CROSS SECTION FOR 12C
C     A. DEUR. JLAB 03/11/03

C
C     ES     - INCIDENT ENERGY (GEV)
C     THETAT - SCATTERING ANGLE (RADIANS)
C     MT     - TARGET MASS (GEV)
C     RECOIL - RECOIL FACTOR
C     EP     - SCATTERED ENERGY (GEV)
C     Q2     - MOMENTUM TRANSFER SQUARED (GEV**2)
C     MOTT   - MOTT CROSS SECTION
C     FC12     - NUCLEAR FORM FACTOR
C
      REAL ES, THETA, MT, RECOIL, EP, Q2, MOTT, FC12
      PARAMETER (MT = 11.178)
C

      RECOIL = 1./(1.+ES/MT*(1.-COS(THETA)))
      EP = RECOIL*ES
      Q2 = 2.*MT*(ES-EP)
      cross_section=RECOIL*MOTT(ES,THETA)*(FC12(Q2)**2.)! in microbarn	
      !PRINT*,es,theta*180./3.1416
      !print*,'Xsec= ',cross_section,' microbarn/sr'
      return
      END
C
C     MOTT CROSS SECTION
C
	REAL FUNCTION MOTT(ES, THETA)
	IMPLICIT NONE
C       
C       ES    - INCIDENT ENERGY (GEV)
C       THETA - SCATTERING ANGLE (RADIANS)
C       Z     - TARGET ATOMIC NUMBER
C       ALPHA - FINE STRUCTURE CONSTANT
C       
	REAL ES, THETA, ALPHA, C2, S,hbc2,Z,ST
	PARAMETER (ALPHA = 1./137.035989561)
	PARAMETER (Z = 6)
	parameter (hbc2=0.38938)
C
	C2 = COS(0.5*THETA)
	ST = SIN(0.5*THETA)
	S = Z*ALPHA*C2/(2.*ES*ST*ST) 
	MOTT = S*S*hbc2*1000.	! microbarn
	RETURN
	END

      real function FC12(QSQP)
C
C     K. Slifer 10/04/02
C     12C Form Factor
C     See K. C. Stansfield et al. PRC 3, 1448 (1971)
C

      HBARC= 197.327053 ! MEV-FM
      ZT   = 6.         ! ATOMIC CHARGE
      AT   = 12.        ! ATOMIC NUMBER
      XALPHA=(ZT-2.)/3.
      QSQP=QSQP*1000000. !now in MeV^2
      Q2FM= QSQP/HBARC**2                   ! fm^-2
      if (Q2FM.LE.3.2) THEN
         xa=1.64                            ! fm
      elseif(Q2FM.GE.3.5) then
         xa=1.68                            ! fm
      else
         xa=0                               ! fm
      endif
      xa=xa/HBARC                           ! 1/MeV

      FC12  = 1.0 - XALPHA/(2.0*(2.+3.*XALPHA)) *QSQP*xa**2
      FC12  = FC12*EXP(-(QSQP*xa**2)/4.)
      if(Q2FM.GE.3.2.and.Q2FM.LE.3.5) then  ! DIFFRACTION MINIMUM
         FC12  = 1.D-5
      endif
      QSQP=QSQP/1000000. !back to GeV^2

      RETURN
      end



****************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function to get the energy loss by external with the right proba due 
!Ionisation and bremstralung. Formula are from Xiaodong Jiang Thesis
! and Claude Marchand Thesis
* CEBAF, A. D 10/25/1999
!Note that the function does not depend much on a and z. Their purpose is 
! just to get a more accurate b.
! input and output Units are GeV and radiation length
      real function incl(E0,z,t)
 
      implicit none
 
      real b,E0,rndm,t,z,rand

      b=4./3.*(1.+(z+1.)/9./(log(1194.*z**(-2./3.))+z*log(
     +184.15*z**(-1./3.))))
      incl=(E0)*(rndm(1)*0.999)**(1./(b*t)) ! Bremstralung

 	if (.not.((incl.ge.-100.).and.(incl.le.10**20.))) then
 	   print*,'Arg, Something wrong with external RC for this event'
           print*,'external:',E0*1000.,z,t,incl
	   incl=10**20.
 	endif		
      end
!routine to compute the probability to lose energy by internal
! Formula are from C Marchand thesis
* CEBAF, A. D 10/25/1999
* Input Units are GeV  and radian
      real function intrad(E0,thp)

      implicit none

      real alpha,E,E0,me,Mc12,nu,pi,rndm,qsq,thp,rand

      alpha=1/137.
      pi=3.141592654
      Mc12=11.178 ! masse nucleus in GeV   
      me=0.000511 ! mass e- in GeV   
 
      E=E0/(1.+2.*E0/Mc12*(sin(thp/2.))**2.)
      qsq=2.*E0*E*(1.-cos(thp))
      nu=2.*alpha/pi*(log(qsq/me**2.)-1.) 
      nu=nu/2.
 	intrad=E*(rndm(1)*0.999)**(1./(nu)) 
	if (.not.((intrad.ge.-100000.).and.(intrad.le.10000000.))) then
	   print*,'Arg, Something wrong with internal RC for this event'
      print*,'internal:',' thp:',thp,' Q^2:',qsq,' rcw:',intrad
	endif	
      end

!function to get energy loss by ionization
!Formula are from the Leo
* CEBAF, A. D 11/27/1999
! input and output Units are GeV and radiation length
      real function ioni(E0,a,z,xd)
 
      implicit none
 
      real a,betas,betasmo,d0,de,E0,ksi
      real lneps,pi,phil,rndm,u,xd,z,rand
      integer j
      parameter (pi=3.1415926)
      common/landau/ld(2,68)
      real ld !landau integral ld(1,i)=lambda, ld(2,i)=int(phi(lambda))
      E0=E0*1000.!convert E0 in MeV
      betas=1-(0.511/E0)**2.
      betasmo=(0.511/E0)**2. ! beta^2 of the electron
      ksi=0.154/betas*z/a*xd ! warning, ksi in MeV
      lneps=log((betasmo)*(0.0135*z)**2./2./0.511/betas)+betas
      d0=ksi*(log(ksi)-lneps+0.37)  ! ionisation energy loss (MeV) most probable      
      u=rndm(1)
      do j=1,68 ! look at the inverse of the integrated landau tail
         if(u.le.ld(2,j)) then ! interpolate
            phil=(ld(1,j)-ld(1,j-1))/(ld(2,j)-ld(2,j-1))
            phil=phil*(u-ld(2,j-1))+ld(1,j-1)
         goto 2
         endif
      enddo
 2    continue
      de=d0+0.06+phil*ksi!express the Landau shape in term of energy loss (MeV)
      E0=E0/1000.
      ioni=de/1000.
     
      end 
!function to get most probable energy loss by ionization
* CEBAF, A. D 11/27/1999
! input and output Units are GeV and radiation length
!      real function deltazero(E0,a,z,xd)

!      implicit none
 
!      real a,betas,betasmo,d0,E0,ksi
!      real lneps,xd,z
 
!      E0=E0*1000.!convert E0 in MeV
!      betas=1-(0.511/E0)**2.
!      betasmo=(0.511/E0)**2. ! beta^2 of the electron
!      ksi=0.154/betas*z/a*xd	! warning, ksi in MeV
!      lneps=log((betasmo)*(0.0135*z)**2./2./0.511/betas)+betas
!      d0=ksi*(log(ksi)-lneps+0.37)  ! ionisation energy loss (MeV) most probable 
!      E0=E0/1000. 	
!      deltazero=d0/1000. !(GeV)
    
!      end 

!! subroutine to fill the common block with landau dist
      subroutine stuffit
	implicit none
	integer i
      common/landau/ld(2,68)
      real lambda(68),intphi(68),ld
      data lambda/-3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-.5,.0 ,.5 ,1.0,1.5 ,2. 
     +,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5 ,9.0 ,9.5,
     +10.0,10.5,11.0 ,11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0,
     +16.5,17.,17.5,18.,18.5,19.,19.5,20.,20.5,21.,21.5,22.,22.5,23.,
     +23.5,24.0,24.5,25.0,25.5,26.,26.5,27.,27.5,28.,28.5,29.,29.5,30. /
      data intphi/.0,5.00000E-04,6.50000E-03,2.90000E-02,8.15000E-02,
     +.1515,.239,.329,.414,.4875,.55,.6025,.6475,.6875001,.7195001,
     +.7455001,.7670001,.7850001,.8005001,.8145001,.8270001,.8385001,
     +.8482501,.8570001,.8650001,.8722501,.8787501,.8847501,.8902501,
     +.8952501,.8998752,.9041252,.9080002,.9115002,.9147502,.9177502,
     +.9205002,.9230002,.9253752,.9276252,.9297502,.9317502,.9336252,
     +.9353752,.9370502,.9386502,.9402002,.9417002,.9430752,.9443253,
     +.9454502,.9464502,.9474002,.9482752,.9491102,.9498602,.9505252,
     +.9511502,.9517627,.9523627,.9529527,.9535328,.9541028,.9546627,
     +.9552028,.9557328,.9562528,.9567627 /
	do i=1,68
	   ld(1,i)=lambda(i)
	   ld(2,i)=intphi(i)
	enddo
      end
!!subroutines for row wise ntuples
      subroutine initpaw
      common/pawc/blanc(150000)
      character*8 tags(25)
*      data tags/'zlab','deltat','thetat','phit','yt','xs','asy'
*     +,'wmm','mott'/
      data tags/'zreact','xs','asy','wmm','mott','yt','phit'
     +,'thetat','deltat','xfoc','yfoc','phfoc','thfoc'
     +,'Eps','Epor','Etemp','ang','x_sieve','y_sieve','angstrag'
     +,'yor','wor','dpor','thor','phor'/
*	data tags/'zlab','xs','asy','wmm','mott','yor','phor'
*     +,'thor','dpp','xfoc','yfoc','phfoc','thfoc'/
      call hlimit (150000)
      call hlimit (150000)
      call hropen (55,'bidon','mchc12ntp.hbook','n',8190,irc)
      call hbookn(1,'truc',25,'bidon',1000,tags)
      end
c
      subroutine endpaw
      call hrout(0,icycle,' ')
      call hrend('bidon')
      end






