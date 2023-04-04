       integer irlen,ilen
       character*9 rankfile
       character*2 sexp
       character*80 comm

       sexp = "04"
       open(10,file='rankfilename',status='old')
       read(10,*) irlen
       read(10,*) rankfile
       close(10)
       comm = "mv "//rankfile(1:irlen)//" "//
     &   rankfile(1:irlen)//"_"//sexp
        write(6,*) comm
c       y=system(comm)
 
        comm = "rm ~/Sci/ECOM/GA/PT/Output/*"
        ilen = 28
        comm = comm(1:ilen)//rankfile(irlen-3:irlen)//"*"
        write(6,*)ilen, comm
 
        end
