      subroutine mksid()
c     converts id array to sid
c     isz = size of id array

      include 'genetic'
      character*1 s0,s1,sidi
      
      s0 = '0'   ! closed cell
      s1 = '1'   ! open cell


      do i=1,nc  ! loop cells

c      determine if open or closed
        if (id(i).eq.0) then 
          sidi=s0 
        else 
          sidi=s1
        endif

c      append string name
        if (i.eq.1) then
          sid = sidi
        else
          sid = sid(1:i-1)//sidi
        endif
      enddo

      write(6,10) (sid)
 10   format(a9)
      return
      end



