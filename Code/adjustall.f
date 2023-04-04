c  adjusts model variables from one time step to another

      subroutine adjustall(dtinew)

      include 'comdeck'
      real adjratio,dtinew

      adjratio = dtinew/dti
c      write(6,*) 'adjusting time step ',dti,dtinew,adjratio
c---------------------------------------------------------      
c   2D prognostic variables
c---------------------------------------------------------      

      call adjustvar2d(ua,uab,adjratio)
      call adjustvar2d(va,vab,adjratio)
      call adjustvar2d(utf,utb,adjratio)
      call adjustvar2d(vtf,vtb,adjratio)
      call adjustvar2d(el,elb,adjratio)
      call adjustvar2d(et,etb,adjratio)
      call adjustvar2d(egf,egb,adjratio)

      write(6,*) 'finished 2d variables'
c---------------------------------------------------------      
c   3D prognostic variables
c---------------------------------------------------------      

c      call adjustvar3d(q,qb,adjratio,fsm)
c      call adjustvar3d(q2f,q2b,adjratio,fsm)
      call adjustvar3d(q2,q2b,adjratio,fsm)
      call adjustvar3d(q2l,q2lb,adjratio,fsm)
      call adjustvar3d(u,ub,adjratio,dum)
      call adjustvar3d(v,vb,adjratio,dvm)
      call adjustvar3d(w,wb,adjratio,fsm)
      call adjustvar3d(s,sb,adjratio,fsm)
      call adjustvar3d(t,tb,adjratio,fsm)
      call adjustvar3d(ccc,cccb,adjratio,fsm)
      
      write(6,*) 'finished 3d variables'
      return
      end
