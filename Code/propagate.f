c  apply fitness testing, survival, mutation and propagation  

       include 'genetic'
c       include 'comdeck'

       nkeep = 2   ! initial # sids to keep in roulette

c      get most recent generation #
       open(10,file='SID_00.dat',status='old')
       read(10,*) sid
       read(10,*) igen0,ttest,tmin
      

       call random_seed()
       call rdfitness(igen0) 
       call roulette(nkeep)
       call mutate2(nkeep)

       end

