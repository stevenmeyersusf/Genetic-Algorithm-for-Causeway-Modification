c      test development of core genetic code

       include 'genetic'
       include 'comdeck'
c       integer nkeep

       nkeep = 2   ! initial # sids to keep in roulette
       igen = 1    ! generation to read and propagate
       ttest = 30000
       tmin = 20000
       call rdfitness(igen) 
       call roulette()
c       write(6,*) 'moving pastd roulette '
       call mutate(nkeep)

       end

