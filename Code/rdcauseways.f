      subroutine rdcauseways(cfile)
c     reads OTB causeway cell coordinates
c     April 12, 2022, S. Meyers

c     nc is total # points
c     idxc, branch index, begins counting in SW corner of OTB
  
      include 'comdeck'
      include 'genetic'

      integer iluni,iycut
      dimension com(80)
      character*45 cfile

      iycut = 51  ! from cutOTB.m 
      ilun = 10
      open(ilun,file=cfile,status='old')
      read(ilun,200) com
      read(ilun,*) nc_all   ! total number of all causeway cells
      if (nc_all.lt.nc) write(6,*) "ERROR in CELL COUNT, rdcauseways.f"

c     initialize ncellspan
      do i=1,nspan
        ncellspan(i) = 0
      enddo
   
c     count cells in each span
      do i=1,nc_all
        read(ilun,*) is(i),js(i),idxc(i)
        ncellspan(idxc(i)) = ncellspan(idxc(i))+1
        js(i) = js(i)-iycut
      enddo
 100  close(ilun)
 200  format(80a1)

      return
      end
      

