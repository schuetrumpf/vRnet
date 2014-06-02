      INTEGER I,N,MA(10000)

      CHARACTER*5 INAM(10000),NAMW(10000),NAMR(6)

      CHARACTER*5 NAME
      INTEGER NA,NB
      DOUBLE PRECISION A,SP,BE,G(24)
      DOUBLE PRECISION ZA(10000),ZP(10000),ZN(10000),BI(10000)
      CHARACTER*72 line

                                  !Chap=chapter number (i.e. rxn type)
      INTEGER Chap,temp,Nnuc(8)   !Total number of nuclide per reaction
      DATA Nnuc/2,3,4,3,4,5,6,5/  !total # = reactants + products

      CHARACTER*4 LABEL
      CHARACTER*1 nr,vw
      DOUBLE PRECISION Q,P(7),I1(6),I2
      
      








      open(1,file='INPUT/sunet',status='old')
      open(2,file='INPUT/winvn',status='old')
      open(3,file='INPUT/rawreaclib',status='old')                          
      open(4,file='netwinv',status='unknown')
      open(5,file='netsu',status='unknown')                                              



C==================================
C     READ IN SUNET
C

      N=0
 11   N=N+1                                                             
      READ(1,99,END=12) INAM(N)                                        
      GOTO 11                                                           
 12   N=N-1                                                             
      WRITE(6,98)(INAM(I),I=1,N)     
                                 
 99   FORMAT(2A5)                                                       
 98   FORMAT(1X,15A5)                                                   

C==================================
C     READ IN WINVN
C

      READ( 2,97) line                                                
      write(4,91) N                                                
      READ( 2,97) line                                                
      write(4,97) line                                                
      NW=0                                                               
      I=1                                                               
 21   NW=NW+1                                                             
      READ(2,96) NAMW(NW)                                              
      IF(NW.EQ.1) GOTO 22                                                
      IF(NAMW(NW).EQ.NAMW(NW-1)) GOTO 23                                  
 22   IF(NAMW(NW).NE.INAM(I)) GOTO 21                                    
      write(4,96) NAMW(NW)                                              
      MA(I)=NW                                                           
      I=I+1                                                             
      GOTO 21                                                            
 23   NW=NW-1                                                             
      I=1                                                               
                                                             
C==================================
C     WRITE FOR NUCLIDES ONLY IN SUNET
C

      DO 30 L=1,N                                                    
 31      K=MA(L)-I+1    
         DO 32 J=1,K
         READ(2,95) NAME,A,NA,NB,SP,BE                                  
 32      READ(2,94) (G(M),M=1,24)                                       
         I=MA(L)+1                                                         
         write(4,95) NAME,A,NA,NB,SP,BE                                  
         write(4,94) (G(M),M=1,24)  
 30   CONTINUE                                                          
  


 97   FORMAT(A72)                                                      
 96   FORMAT(2A5)                                                       
 95   FORMAT(A5,F12.3,2I4,F6.1,F10.3)                                   
 94   FORMAT(8F9.2)               
      
C==================================
C     
C



C********************************************************************
C     READING FFN WEAK RATES WOULD GO HERE 
C     SO WE CAN REPLACE SIMPLE BETA-DECAYS
C     WITH TEMP & DENS DEPENDENT RATES
C     FOR NOW, WE IGNORE

      NM=0                                                              
C
C********************************************************************
      nline=0
      WRITE(5,91) NM,NM                                                  

 91   FORMAT(2I5)

 40   READ(3,93,END=41,ERR=42) K,(NAMR(J),J=1,6),LABEL,nr,vw,Q                   
      READ(3,92,ERR=43) (P(J),J=1,7)                                          
      nline=nline+3
      if(K.NE.0) then
         Chap=K     
         write(5,93)K,(NAMR(J),J=1,6),LABEL,nr,vw,Q    
         WRITE(5,92)(P(J),J=1,7)                                          
         goto 40
      endif
      goto 44
 42   write(6,*) 'ERROR: header line', nline
      stop
 43   write(6,*) 'ERROR: body line', nline
      write(6,93) K,(NAMR(J),J=1,6),LABEL,nr,vw,Q                   
      write(6,92) (P(J),J=1,7)                                          
      stop
 44   continue
 
      do J=1,Nnuc(Chap)
         I1(J)=0.
      enddo

      do I=1,N
      do J=1,Nnuc(Chap)
         if(adjustl(NAMR(J)).eq.adjustl(INAM(I)))I1(J)=I1(J)+1.
      enddo
      enddo

      temp=Nnuc(Chap)
      if((Chap.eq.8).and.(adjustl(NAMR(5)).eq."     "))temp=4
      I2=1.
      do J=1,temp
         I2=I2*I1(J)
      enddo

      if(I2.eq.1.)then
         WRITE(5,93) K,(NAMR(J),J=1,6),LABEL,nr,vw,Q                         
         WRITE(5,92)(P(J),J=1,7)                                          
      endif

 93   FORMAT(I1,4X,6A5,8X,A4,a1,a1,3X,1pE12.5)                          
 92   FORMAT(4E13.6)                                                    

      goto 40
 41   continue



      end
