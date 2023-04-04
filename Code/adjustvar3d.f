c     adjusts 3d variables to new time step
      subroutine adjustvar3d(var,varb,adjratio)

      real adjratio,var(im,jm,kb),varb(im,jm,kb)

         do 100 i=1,im
         do 200 j=1,jm
            do 300 k=1,kb
               dvij = (var(i,j,k)-varb(i,j,k))*adjratio
               varb(i,j,k) = var(i,j,k)-dvij
 300        continue
 200     continue
 100     continue

      return
      end
