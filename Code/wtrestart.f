      subroutine wtrestart(ADVVA,ADVUA, ADVVV, ADVUU)

      include 'comdeck'

         WRITE(IUWRS) 
     .     intcold,INT,DZR,Z,ZZ,DZ,DZZ,H,H1,H2,D,DT,ANG,AAM,AAM2D,
     .     ART,ARU,ARV,DUM,DVM,FSM,COR,CURV42D,WUBOT,WVBOT,
     .     UA,UAB,VA,VAB,EL,ELB,ETF,ET,ETB,EGF,EGB,UTF,UTB,
     .     VTF,VTB,ADVUU,ADVVV,ADVUA,ADVVA,KM,KH,KQ,Q2,Q2B,
     .     Q2L,Q2LB,L,U,UB,W,V,VB,T,TB,S,SB,RHO,RMEAN,TMEAN,SMEAN,
     .     CCC,CCCB,CCCMEAN
         
         CLOSE (IUWRS)

         return


         end
