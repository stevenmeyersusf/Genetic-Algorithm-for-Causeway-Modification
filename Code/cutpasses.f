       subroutine cutpasses()
c      changes model depth matrix h

       include 'genetic'
       include 'comdeck'
       integer ic
       integer*4 id
       character*1 sidi
       real hinterp

c-------------------------------------------------
c  adjust H and U,V masks
c-------------------------------------------------
c       write(6,*) 'adapting causeways...'
       do i=1,nc
         sidi(1:1) = sid(i:i)
c         write(6,*) i,'s=',sidi,is(i),js(i)
         if (sidi.eq.'1') then !  is a cut
             ii = is(i)
             jj = js(i)
             hinterp = (h(ii,jj-1)+h(ii,jj+1))/2  ! get avg H above and below
             h(ii,jj) = hinterp
c             fsm(ii,jj) = 1  ! now in setdom.f
c             dum(ii,jj) = 1
c             dvm(ii,jj) = 1
c             write(6,*) ii,jj,h(ii,jj-1),h(ii,jj+1),h(ii,jj)
         endif
       enddo

c-------------------------------------------------------
c   check U mask and adjust if needed
c   only checking U mask bc causeway cells are 
c   single-file along-x
c-------------------------------------------------------
c      do i=1,nc
c         sidi(1:1) = sid(i:i)
c         if (sidi.eq.'1') then !  is a cut
c           ii = is(i)
c           jj = js(i)  
c           IF(FSM(II,JJ).EQ.0.0.AND.FSM(II+1,JJ).NE.0.0) DUM(II+1,JJ)=0.0
c         endif
c       enddo


       return
       end
         
       


