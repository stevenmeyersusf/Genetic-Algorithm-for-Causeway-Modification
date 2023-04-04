      subroutine testcfl(cflmin)
c     tests basic CFL condition across grid

      include 'comdeck'
      real  umax,vmax,utime,vtime,r1,r2
      integer imin,jmin

c      write(6,*) 'testing CFL ',int
      cflmin = 99999999.
      do 100 i=1,im
      do 100 j=1,jm
         if (fsm(i,j).eq.0) goto 100  ! skip land points
         umax = 0.
         vmax = 0
         do 110 k=1,kb
            uk = u(i,j,k)*dum(i,j)
            if ((abs(uk).gt.abs(umax)).and.(uk.gt.0)) umax=abs(uk)
            uk = u(i+1,j,k)*dum(i+1,j)
            if ((abs(uk).gt.abs(umax)).and.(uk.lt.0)) umax=abs(uk)

            vk = v(i,j,k)*dvm(i,j)
            if ((abs(vk).gt.abs(vmax)).and.(vk.gt.0)) vmax=abs(vk)
            vk = v(i+1,j,k)*dvm(i+1,j)
            if ((abs(vk).gt.abs(vmax)).and.(vk.lt.0)) vmax=abs(vk)
 110     continue
         utime = h1(i,j)/umax
         vtime = h2(i,j)/vmax
         if ((utime.lt.vtime).and.(utime/dti0.lt.cflmin)) then 
            cflmin=utime/dti0
            imin=i
            jmin=j
         endif
         if ((vtime.le.utime).and.(vtime/dti0.lt.cflmin)) then 
            cflmin=vtime/dti0
            imin=i
            jmin=j
         endif
 100  continue

      if (cflmin.lt.1.5) write(6,*)
     &     'CFL WARNING ',int,int/60.,imin,jmin,cflmin
      
c      if (mod(int,60*24.eq.0)
c      write(6,*) 'cflmin MIN CFL ',imin,jmin,cflmin,el(imin,jmin),wsp

c      write(6,*) 'finished testing CFL '

      return
      end
