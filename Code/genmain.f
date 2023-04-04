c23456789012345678901234567890123456789012345678901234567890123456789012

c  trial main routine for genetic development
c  when finished, must manually rename rank* to rank*_exp#
c       routine main

      include 'comdeck'
      include 'genetic'

      character*1 a1
      character*2 sexp
      character*45 gfile
      character*80 comm
      integer iexp,ilen 

      call random_seed()

c     must be consistent with definition in main.f
      data ispan/0,0,0,0,1,0/
      data ncutx/0,0,0,0,2,0/
      data is/8,9,10,11,12,13,14,15,16/  ! define span cells
      data js/33,32,32,32,32,32,32,32,32/
      data idxc/5,5,5,5,5,5,5,5,5/
      data ncellspan/0,0,0,0,9,0/
      data ttest/70560./   ! INT time to do fitness test (ftest0.f)
      data tmin/30240./     ! INT time to calc n0 in ftest0.f

c#################################################
c    BEGIN Experiment Loop
c#################################################
      do iexp=1,nexp
        if (iexp.lt.10) then 
           write(a1,'(i1)') iexp
           sexp = "0"//a1
        endif
        if (iexp.ge.10) then 
           write(sexp,'(i2)') iexp
        endif
c---------------------------------------------
c 1. make random grids for generation 1
c 2. run GA: bash+FORTRAN code
c---------------------------------------------
        x = system('rm SID*dat')  ! remove priors to avoid confusion
        npop = 4  ! initialize # active in population
        igen = 1  ! first generation
        write(6,*) 'genmain igen npop',igen,npop
        gfile = '/home/meyers/Sci/ECOM/GA/PT/causeway_grid.txt'
        call mkrandgrid(igen) ! one-time call, make initial population

c       begin GA 
        x=system('sh modelgroup')
 
c--------------------------------------------------------
c     rename rank file with number of experiment
c--------------------------------------------------------
c       get rank file name created by 'modelgroup'
        open(10,file='rankfilename',status='old')
        read(10,*) irlen
        read(10,*) rankfile
        close(10)
        comm = "mv "//rankfile(1:irlen)//" "//
     &    rankfile(1:irlen)//"_"//sexp 
        y=system(comm)

c--------------------------------------------------------
c    delete prior rank file for this generation
c--------------------------------------------------------
        
c        comm = "rm ~/Sci/ECOM/GA/PT/Output/*"
c        ilen = 28
c        comm = comm(1:ilen)//rankfile(irlen-4:irlen)//"*"
c        z = system(comm)
      enddo   ! loop experiments 


      end
