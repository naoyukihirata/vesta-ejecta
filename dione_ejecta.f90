!This program calculates the landing locations of the ejecta particles launched from the Creusa crater of Dione.  
!This program is made by Naoyuki Hirata, Kobe Univ. 
!This program follows the method described by N. Hirata et al. [2021], Icarus 354, 114073.
!The initial launch velocity of a particle is given by Housen and Holsapple [2011], Icarus 211, 856–875.
!The equation of motion is performed in a rotating reference frame (i.e. the asteroid-fixed frame) described by Scheeres et al. [2002], Asteroid III, Univ. of Arizona Press. 

     implicit none
     integer :: i, k,j
     integer, parameter :: num = 10000
     real(8) :: c_lon, c_lat,pi,rtem1,rtem2
     integer :: loop,azimuth,Target, Frame
     real(8) :: pp(3),pv(3),dx(3),dp(3),dx_prev(3),dp_prev(3),dt,g,launch_point
     real(8) :: gm,distance,w,w2, radius,launch,rotation,c_radius
     real(8) ::  c_nz,c_nx,c_ny,Q,c_Q,c_cQ,c_sQ,im_lat,im_lon
     real(8) :: mcx,mcy,mcz,cy,cx,cz,cQ,sQ,nz,ny,nx,polar_lon,polar_lat
     real(8) :: c1, h, u, n2, kk, ro, yy, p
     real(8) :: mass_asteroid,  mass_sun, semi_major_axis,hill_radius, angle

   !!!input parameter

   pi = acos(-1.0d0)
   angle = 30.                           !initial launch angle of ejecta particle from horizontal
   azimuth = 180.                        !initial launch azimuth (from 0 to 360 degree)  
   dt = 1.d0                             !time step
   radius = 561.4d0*1000.d0              !radius of the asteroid (Dione) [m]
   mass_asteroid   =1.095d21             !mass of the asteroid (Dione) [kg]
   mass_sun = 5.688d26                   !mass of the center object (Saturn) [kg]
   semi_major_axis = 377400.*1000.d0     !s.m.a. of the asteroid [m]
   !rotation = 5.342                      !rotation period of the asteroid [h]
   c_radius = 0.5d0*dble(36.2*1000.)     ! crater radius (Creusa crater) [m]
   c_lat = dble(49) *pi/180.d0       !latitude of the crater center
   c_lon = dble(284) *pi/180.d0          !east longitude of the crater center

   open(11,file="result.txt")  ! calculation result output file
    write(11,*) "# Vej, x/n2Rc, lat. of landing location, east long., travelling time"

   Frame = 2
   ! 0 : non-rotating reference frame
   ! 1 : rotating reference frame (asteroid)
   ! 2 : hill's coordinate system (synchronously-rotating satellite)
   ! 3 : non-synchronously-rotating satellite. In this case, rotation = (synodic rotation period ) in Line 29

   Target = 4       ! Target type is selected from 1 to 8. Parameters are from Housen and Holsapple [2011]

   if (Target .eq. 1) then
   u  = 0.55
   kk = 0.2
   c1 = 1.5
   h  = 0.68
   n2 = 1.5
   p  = 0.5
   ro = 1000.
   endif
   if (Target .eq. 2) then
   u  = 0.55
   kk = 0.3
   c1 = 1.5
   h  = 1.1
   n2 = 1
   p  = 0.5
   ro = 3000.
   yy = 30
   endif
   if (Target .eq. 3) then
   u  = 0.46
   kk = 0.3
   c1 = 0.18
   h  = 0.38
   n2 = 1
   p  = 0.3
   ro = 2600.
   yy = 0.45
   endif
   if ((Target .eq. 4) .or. (Target .eq. 5)) then
   u  = 0.41
   kk = 0.3
   c1 = 0.55
   h  = 0.59
   n2 = 1.3
   p  = 0.3
   ro = 1600.
   endif
   if (Target .eq. 6) then
   u  = 0.45
   kk = 0.5
   c1 = 1.
   h  = 0.8
   n2 = 1.3
   p  = 0.3
   ro = 1500.
   endif
   if (Target .eq. 7) then
   u  = 0.4
   kk = 0.3
   c1 = 0.55
   h  = 0.4
   n2 = 1.
   p  = 0.3
   ro = 1500.
   yy = 4.e-3
   endif
   if (Target .eq. 8) then
   u  = 0.35
   kk = 0.32
   c1 = 0.6
   h  = 0.81
   n2 = 1.
   p  = 0.2
   ro = 1200.
   yy = 2.e-3
   endif


      gm = mass_asteroid*6.67408d-11  
      g = gm/(radius*radius)
      hill_radius = semi_major_axis * (mass_asteroid/(mass_sun*3.d0))**(1.d0/3.d0)
      k = 0  
      c_radius = c_radius/n2
      Q =  dble(azimuth)*2.d0*pi/360
      cQ = dcos(Q)
      sQ = dsin(Q)

      if (Frame .eq. 0) then
        w = 0.d0                              !rad/s
      elseif (Frame .eq. 1) then
        w = 2.*pi/(rotation*3600.)            !rad/s
      elseif ((Frame .eq. 2) .or.(Frame .eq. 3)) then
        w = dsqrt(mass_sun*6.67408d-11/(semi_major_axis**3))
        w2 = 2.*pi/(rotation*3600.)
        c_lon=c_lon-pi
      endif

      polar_lon = 0.d0*pi/180.d0
      polar_lat = 0.d0*pi/180.d0
      nz = dsin(polar_lat)
      nx = dcos(polar_lon)*dcos(polar_lat)
      ny = dsin(polar_lon)*dcos(polar_lat)


   do loop = 1, num-1

    !initial launch position and velocity of ejecta particles

    launch_point = n2*c_radius*dble(loop)/(radius*dble(num)) 
    if ((Target .eq. 1) .or.(Target .eq. 4) .or.(Target .eq. 5) .or.(Target .eq. 6)) then 
      launch = c1*(h*(4.*pi/3.)**(1./3.))**((-2.-u)/(2.*u))*dsqrt(g*c_radius) &
        *(launch_point/(c_radius/radius))**(-1./u) &
        *(1.d0 -  launch_point/(n2*c_radius/radius) )**p
    else
      launch = c1*(h*(4.*pi/3.)**(1./3.))**(-1./u)*dsqrt(yy/ro) &
        *(launch_point/(c_radius/radius))**(-1./u) &
        *(1.d0 - launch_point/(n2*c_radius/radius) )**p
    endif
    pp(1) = radius*dcos(launch_point)
    pp(2) = radius*dsin(launch_point)
    pp(3) = 0.d0
     rtem2 = atan2(pp(2), pp(1))
    pv(1) = launch *cos(rtem2+0.5*pi-angle*pi/180.)
    pv(2) = launch *sin(rtem2+0.5*pi-angle*pi/180.)
    pv(3) = 0.d0 

    !using representation matrix of Rodrigues' rotation formula, initial launch position and velocity were obtained

    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( cQ+nx*nx*(1.d0-cQ) )*cx      &
           +(nx*ny*(1.d0-cQ)-nz*sQ)*cy      &
           +(nx*nz*(1.d0-cQ)+ny*sQ)*cz
      mcy = (ny*nx*(1.d0-cQ)+nz*sQ)*cx      &
           +(cQ+ny*ny*(1.d0-cQ)  )*cy      &
           +(ny*nz*(1.d0-cQ)-nx*sQ)*cz
      mcz = (nz*nx*(1.d0-cQ)-ny*sQ)*cx      &
           +(nz*ny*(1.d0-cQ)+nx*sQ)*cy      &
           +(cQ+nz*nz*(1.d0-cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz

    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( cQ+nx*nx*(1.d0-cQ) )*cx      &
           +(nx*ny*(1.d0-cQ)-nz*sQ)*cy      &
           +(nx*nz*(1.d0-cQ)+ny*sQ)*cz
      mcy = (ny*nx*(1.d0-cQ)+nz*sQ)*cx      &
           +(cQ+ny*ny*(1.d0-cQ)  )*cy      &
           +(ny*nz*(1.d0-cQ)-nx*sQ)*cz
      mcz = (nz*nx*(1.d0-cQ)-ny*sQ)*cx      &
           +(nz*ny*(1.d0-cQ)+nx*sQ)*cy      &
           +(cQ+nz*nz*(1.d0-cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz

      c_nz = 0.d0 
      c_nx = 0.d0 
      c_ny = 1.d0 
      c_Q = -c_lat
      c_cQ = dcos(c_Q)
      c_sQ = dsin(c_Q)
    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz
    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz
      c_nz = 1.d0
      c_nx = 0.d0 
      c_ny = 0.d0 
      c_Q = c_lon
      c_cQ = dcos(c_Q)
      c_sQ = dsin(c_Q)
    cx = pp(1)
    cy = pp(2)
    cz = pp(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pp(1) = mcx
    pp(2) = mcy
    pp(3) = mcz
    cx = pv(1)
    cy = pv(2)
    cz = pv(3)
      mcx = ( c_cQ+c_nx*c_nx*(1.d0-c_cQ) )*cx      &
           +(c_nx*c_ny*(1.d0-c_cQ)-c_nz*c_sQ)*cy      &
           +(c_nx*c_nz*(1.d0-c_cQ)+c_ny*c_sQ)*cz
      mcy = (c_ny*c_nx*(1.d0-c_cQ)+c_nz*c_sQ)*cx      &
           +(c_cQ+c_ny*c_ny*(1.d0-c_cQ)  )*cy      &
           +(c_ny*c_nz*(1.d0-c_cQ)-c_nx*c_sQ)*cz
      mcz = (c_nz*c_nx*(1.d0-c_cQ)-c_ny*c_sQ)*cx      &
           +(c_nz*c_ny*(1.d0-c_cQ)+c_nx*c_sQ)*cy      &
           +(c_cQ+c_nz*c_nz*(1.d0-c_cQ) )*cz
    pv(1) = mcx
    pv(2) = mcy
    pv(3) = mcz

   if (Frame .eq. 3) then
    pv(1) = pv(1) -w2*pp(2)
    pv(2) = pv(2) +w2*pp(1)
   endif

 
    !time development using Adams-Bashforth Method

    distance =sqrt( pp(1)**2.+pp(2)**2.+pp(3)**2.)
    rtem1 = dt*gm/distance**3.
    dp(1) =   -pp(1)*rtem1
    dp(2) =   -pp(2)*rtem1
    dp(3) =   -pp(3)*rtem1
    if (Frame .eq. 1) then
      dp(1) = dp(1) +dt*2.*w*pv(2) +dt*w*w*pp(1)
      dp(2) = dp(2) -dt*2.*w*pv(1) +dt*w*w*pp(2)
     elseif ((Frame .eq. 2) .or.  (Frame .eq. 3)) then
      dp(1) = dp(1) +dt*2.*w*pv(2)+dt*3.d0*w*w*pp(1)
      dp(2) = dp(2) -dt*2.*w*pv(1) 
      dp(3) = dp(3) -dt*w*w*pp(3)
     endif
    dx(1) = pv(1)*dt
    dx(2) = pv(2)*dt
    dx(3) = pv(3)*dt

    do i = 1, 3
    pv(i)=pv(i)+dp(i)
    pp(i)=pp(i)+dx(i)
    dp_prev(i) = dp(i)
    dx_prev(i) = dx(i)
    enddo

    do j = 1, int(1000000000./dt)
       distance =sqrt( pp(1)**2.+pp(2)**2.+pp(3)**2.)
       if (distance .gt.  hill_radius ) exit
       rtem1 = dt*gm/distance**3.
       dp(1) =   -pp(1)*rtem1 
       dp(2) =   -pp(2)*rtem1
       dp(3) =   -pp(3)*rtem1
    if (Frame .eq. 1) then
      dp(1) = dp(1) +dt*2.*w*pv(2) +dt*w*w*pp(1)
      dp(2) = dp(2) -dt*2.*w*pv(1) +dt*w*w*pp(2)
     elseif ((Frame .eq. 2) .or.  (Frame .eq. 3)) then
      dp(1) = dp(1) +dt*2.*w*pv(2)+dt*3.d0*w*w*pp(1)
      dp(2) = dp(2) -dt*2.*w*pv(1) 
      dp(3) = dp(3) -dt*w*w*pp(3)
     endif
       dx(1) = pv(1)*dt
       dx(2) = pv(2)*dt
       dx(3) = pv(3)*dt
         do i = 1, 3
         pv(i)=pv(i)+1.5*dp(i)-0.5*dp_prev(i)
         pp(i)=pp(i)+1.5*dx(i)-0.5*dx_prev(i)
         dp_prev(i) = dp(i)
         dx_prev(i) = dx(i)
         enddo

       !collision judgement 
       if (distance .le. radius ) then
         im_lat = atan(pp(3) / sqrt(pp(1)**2.d0+pp(2)**2.d0 ) )*180.d0/pi
         im_lon = atan2(pp(2), pp(1)) *180.d0/pi
            if (Frame .eq. 2) then
            im_lon=im_lon-180.
            elseif (Frame .eq. 3) then
            im_lon=im_lon-180.
            im_lon=im_lon-(dble(j)*dt*w2)*180./pi
            endif
            do
             if (im_lon .le. 0) then
              im_lon=im_lon+360
             else
              exit
             endif
            enddo
         write(11,*) launch,launch_point, im_lat,im_lon,j*dt
         exit
       endif

    enddo
  enddo 

  close(11)
  end






