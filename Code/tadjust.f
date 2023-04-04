c adjusts time step to adapt to CFL stability
      subroutine tadjust(nadj,dint,cflmin,cflmin0)
      include 'comdeck'
      integer idir,nadjnew,i0,i1,i2
      real dtinew,dint0,cavg,thresh(5),dts(6),cfl2

c      data thresh/0.5, 1.00, 6.0,10., 20./     ! cfl thresholds
      data thresh/0.5, 1.00, 4.0, 8., 20./     ! cfl thresholds
      data dts/0.1250, 0.25, 0.5, 1., 2.0, 4./ ! dint steps

      dint0=dint
      cfl2 = (cflmin+cflmin0)/2.

c     find current index
      do 100 i=1,6
         if (dint0.eq.dts(i)) i0=i
 100  continue

c     find next index, must be at least X% below threshold
      i1 = 6
      i2 = 6
      do 200 i=i1-1,1,-1
c         if (cfl2 .lt.thresh(i)) i1=i
         if ((cflmin-thresh(i))/thresh(i).lt.0.01) i1=i
         if ((cflmin0-thresh(i))/thresh(i).lt.0.01) i2=i
 200  continue

c      if (int.gt.4000)write(6,*)cflmin,cflmin0,i1,i2

c     must have strong change in cflmin to trigger change in time step
c      if (abs(cflmin-cflmin0)/cflmin0.lt.0.05) goto 1000  

      if (i1.ne.i2) goto 1000

c     only change one time step level
c      if (i1.gt.i0) i1 = i0+1
c      write(6,*) 'new i ',i0,i1,dint0,cflmin
c      if (i1.lt.i0) i1 = i0-1
c      write(6,*) 'new i ',i0,i1,dint0,cflmin

c    select new 
      if (i1.eq.1) then
         dint = dts(1)
         dtinew = dint*dti0
         nadjnew = 8
      endif

      if (i1.eq.2) then 
         dint = dts(2)
         dtinew = dint*dti0
         nadjnew = 4
      endif         

      if (i1.eq.3) then 
         dint = dts(3)
         dtinew = dint*dti0
         nadjnew = 2
      endif

      if (i1.eq.4) then 
         dint = dts(4)
         dtinew = dint*dti0
         nadjnew = 1
      endif

      if (i1.eq.5) then 
         dint = dts(5)
         dtinew = dint*dti0
         nadjnew = 1
      endif

      if (i1.eq.6) then 
         dint = dts(6)
         dtinew = dint*dti0
         nadjnew = 1
      endif
      
c      write(6,*) '2dti ',dti,dtinew,int
      if (dti.ne.dtinew) then 
c         write(6,*) 'CFL ADJ ',int,thour,dti,dtinew,i0,i1,cflmin,cflmin0
         call adjustall(dtinew)
c         cflmin0 = cflmin
         dti = dtinew
         nadj = nadjnew
      endif

 1000 return
      end
      
