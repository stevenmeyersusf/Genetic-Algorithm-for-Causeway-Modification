c     adjusts 2d variables to new time step
      subroutine adjustvar2d(var,varb,adjratio)

      include 'comdeck'
      real adjratio,var(im,jm),varb(im,jm),dvij

      do 100 i=1,im
         do 200 j=1,jm
            dvij = (var(i,j)-varb(i,j))*adjratio
            varb(i,j) = var(i,j)-dvij
 200     continue
 100     continue

      return
      end
            
