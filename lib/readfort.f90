module readfort

    use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_PTR, C_LOC, c_float, c_long
    
    ! implicit none
    
    contains
    
    function c_to_f_string(s) result(str)
      use iso_c_binding
      character(kind=c_char,len=1), intent(in) :: s(*)
      character(len=:), allocatable :: str
      integer i, nchars
      i = 1
      do
         if (s(i) == c_null_char) exit
         i = i + 1
      end do
      nchars = i - 1  ! Exclude null character from Fortran string
      allocate(character(len=nchars) :: str)
      str = transfer(s(1:nchars), str)
    end function c_to_f_string
    
    !** Convert a Fortran string to a C string
    function f_to_c_string(f_string) result(c_string)
        character(len=*), intent(in) :: f_string
        character(len=1, kind=c_char) :: c_string(len_trim(f_string)+1)
        type(C_PTR) :: str
        integer :: N, i
    
        N = len_trim(f_string)
        do i = 1, N
            c_string(i) = f_string(i:i)
        end do
        c_string(n + 1) = c_null_char
        str = transfer(c_string,str)
    end function f_to_c_string
    

    subroutine readfile12(&
        s, NLINESi,NBUFFi,CGFi,NELIONi,ELOi,GAMMARi,GAMMASi,GAMMAWi) bind(c, name='readfile12')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(in), value :: NLINESi

        integer(c_int), intent(out) :: NBUFFi(NLINESi)
        real(c_float), intent(out) :: CGFi(NLINESi)
        integer(c_int), intent(out) :: NELIONi(NLINESi)
        real(c_float), intent(out) :: ELOi(NLINESi)
        real(c_float), intent(out) :: GAMMARi(NLINESi)
        real(c_float), intent(out) :: GAMMASi(NLINESi)
        real(c_float), intent(out) :: GAMMAWi(NLINESi)

        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)
    
        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        I = 1
        do while (1.eq.1)
            read(1,end=100) NBUFF,CGF,NELION,ELO,GAMMAR,GAMMAS,GAMMAW
            ! IF(I.GE.NLINESi-100) WRITE(6,*)  NBUFF,CGF,NELION,ELO,GAMMAR,GAMMAS,GAMMAW
            NBUFFi(I) = NBUFF
            CGFi(I) = CGF
            NELIONi(I) = NELION
            ELOi(I) = ELO
            GAMMARi(I) = GAMMAR
            GAMMASi(I) = GAMMAS
            GAMMAWi(I) = GAMMAW
            I = I + 1
        end do
100     continue
        close(unit=1)
    end subroutine readfile12

    subroutine writefile12(&
        s, NLINESi,NBUFFi,CGFi,NELIONi,ELOi,GAMMARi,GAMMASi,GAMMAWi) bind(c, name='writefile12')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(in), value :: NLINESi

        integer(c_int), intent(in) :: NBUFFi(NLINESi)
        real(c_float),  intent(in) :: CGFi(NLINESi)
        integer(c_int), intent(in) :: NELIONi(NLINESi)
        real(c_float),  intent(in) :: ELOi(NLINESi)
        real(c_float),  intent(in) :: GAMMARi(NLINESi)
        real(c_float),  intent(in) :: GAMMASi(NLINESi)
        real(c_float),  intent(in) :: GAMMAWi(NLINESi)

        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)
        INTEGER I

        open(UNIT=1,FILE=c_to_f_string(s),STATUS='NEW',FORM='UNFORMATTED')
        do I=1,NLINESi
            NBUFF = NBUFFi(I)
            CGF = CGFi(I)
            NELION = NELIONi(I)
            ELO = ELOi(I)
            GAMMAR = GAMMARi(I)
            GAMMAS = GAMMASi(I)
            GAMMAW = GAMMAWi(I)
            ! IF(I.GE.NLINESi-100) WRITE(6,*)  NBUFF,CGF,NELION,ELO,GAMMAR,GAMMAS,GAMMAW
            write(1) NBUFF,CGF,NELION,ELO,GAMMAR,GAMMAS,GAMMAW
        end do
        close(unit=1)
    end subroutine writefile12

    subroutine readfile19(&
        s, NLINESi,WLVACi,ELOi,GFi,NBLOi,NBUPi,NELIONi,TYPEi,NCONi,NELIONXi,&
        GAMMARi,GAMMASi,GAMMAWi,NBUFFi,LIMi) bind(c, name='readfile19')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(in), value :: NLINESi

        real(c_double), intent(out) :: WLVACi(NLINESi)
        real(c_float), intent(out) :: ELOi(NLINESi)
        real(c_float), intent(out) :: GFi(NLINESi)
        integer(c_int), intent(out) :: NBLOi(NLINESi)
        integer(c_int), intent(out) :: NBUPi(NLINESi)
        integer(c_int), intent(out) :: NELIONi(NLINESi)
        integer(c_int), intent(out) :: TYPEi(NLINESi)
        integer(c_int), intent(out) :: NCONi(NLINESi)
        integer(c_int), intent(out) :: NELIONXi(NLINESi)
        real(c_float), intent(out) :: GAMMARi(NLINESi)
        real(c_float), intent(out) :: GAMMASi(NLINESi)
        real(c_float), intent(out) :: GAMMAWi(NLINESi)
        integer(c_int), intent(out) :: NBUFFi(NLINESi)
        integer(c_int), intent(out) :: LIMi(NLINESi)

        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        INTEGER TYPE
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)
        
        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        I = 1
        do while (1.eq.1)
            read(1,end=100) WLVAC,ELO,GF,NBLO,NBUP,NELION,TYPE,NCON,NELIONX,GAMMAR,GAMMAS,GAMMAW,NBUFF,LIM
            ! WRITE(6,*) 'R19',I,WLVAC,ELO,GF,NBLO,NBUP,NELION,TYPE,NCON,NELIONX,GAMMAR,GAMMAS,GAMMAW,NBUFF,LIM
            WLVACi(I) = WLVAC
            ELOi(I) = ELO
            GFi(I) = GF
            NBLOi(I) = NBLO
            NBUPi(I) = NBUP
            NELIONi(I) = NELION
            TYPEi(I) = INT(TYPE)
            NCONi(I) = NCON
            NELIONXi(I) = NELIONX
            GAMMARi(I) = GAMMAR
            GAMMASi(I) = GAMMAS
            GAMMAWi(I) = GAMMAW
            NBUFFi(I) = NBUFF
            LIMi(I) = LIM
            I = I + 1
        end do
100     continue
 

        close(unit=1)
    end subroutine readfile19

    subroutine writefile19(&
        s, NLINESi,WLVACi,ELOi,GFi,NBLOi,NBUPi,NELIONi,TYPEi,NCONi,NELIONXi,&
        GAMMARi,GAMMASi,GAMMAWi,NBUFFi,LIMi) bind(c, name='writefile19')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(in), value :: NLINESi

        real(c_double),  intent(in) :: WLVACi(NLINESi)
        real(c_float),  intent(in) :: ELOi(NLINESi)
        real(c_float),  intent(in) :: GFi(NLINESi)
        integer(c_int), intent(in) :: NBLOi(NLINESi)
        integer(c_int), intent(in) :: NBUPi(NLINESi)
        integer(c_int), intent(in) :: NELIONi(NLINESi)
        integer(c_int), intent(in) :: TYPEi(NLINESi)
        integer(c_int), intent(in) :: NCONi(NLINESi)
        integer(c_int), intent(in) :: NELIONXi(NLINESi)
        real(c_float),  intent(in) :: GAMMARi(NLINESi)
        real(c_float),  intent(in) :: GAMMASi(NLINESi)
        real(c_float),  intent(in) :: GAMMAWi(NLINESi)
        integer(c_int), intent(in) :: NBUFFi(NLINESi)
        integer(c_int), intent(in) :: LIMi(NLINESi)

        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        INTEGER TYPE
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)

        open(UNIT=1,FILE=c_to_f_string(s),STATUS='NEW',FORM='UNFORMATTED')
        do I=1,NLINESi
            WLVAC = WLVACi(I)
            ELO = ELOi(I)
            GF = GFi(I)
            NBLO = NBLOi(I)
            NBUP = NBUPi(I)
            NELION = NELIONi(I)
            TYPE = TYPEi(I)
            NCON = NCONi(I)
            NELIONX = NELIONXi(I)
            GAMMAR = GAMMARi(I)
            GAMMAS = GAMMASi(I)
            GAMMAW = GAMMAWi(I)
            NBUFF = NBUFFi(I)
            LIM = LIMi(I)
            ! WRITE(6,*) 'W19',I,WLVAC,ELO,GF,NBLO,NBUP,NELION,TYPE,NCON,NELIONX,GAMMAR,GAMMAS,GAMMAW,NBUFF,LIM
            write(1) WLVAC,ELO,GF,NBLO,NBUP,NELION,TYPE,NCON,NELIONX,GAMMAR,GAMMAS,GAMMAW,NBUFF,LIM
        end do
 
        close(unit=1)
    end subroutine writefile19

    subroutine readfile14_20(&
        s, NLINESi,&
        WLi,Ei,EPi,&
        LABELi,LABELx,LABELPi,LABELPx,&
        ISHIFTi,ISHIFTPi,IXFIXFPi,LINESIZEi,AUTOi,OTHER1x,OTHER2x,&
        WLVACi,CENTERi,CONCENi,NELIONi,GAMMARi,GAMMASi,GAMMAWi,REFi,REFx,&
        NBLOi,NBUPi,ISO1i,X1i,ISO2i,X2i,GFLOGi,XJi,XJPi,CODEi,ELOi,GFi,GSi,GRi,GWi,&
        DWLi,DGFLOGi,DGAMMARi,DGAMMASi,DGAMMAWi,DWLISOi,ISOSHIFTi,EXTRA3i) bind(c, name='readfile14_20')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)
        character(len=25) :: SLABEL
      
        integer(c_int), intent(in), value :: NLINESi
            
        real(c_double), intent(out) :: WLi(NLINESi)
        real(c_double), intent(out) :: Ei(NLINESi)
        real(c_double), intent(out) :: EPi(NLINESi)
        character(kind=c_char,len=1),  intent(out) :: LABELi(10,NLINESi)
        character(kind=c_char,len=1),   intent(out) :: LABELPi(10,NLINESi)
        real(c_double), intent(out) :: LABELx(2,NLINESi)
        real(c_double), intent(out) :: LABELPx(2,NLINESi)
        ! character(kind=c_char,len=1),   intent(out) :: OTHER1i(10,NLINESi)
        ! character(kind=c_char,len=1),   intent(out) :: OTHER2i(10,NLINESi)

        integer(c_int), intent(out)  :: ISHIFTi(NLINESi)
        integer(c_int), intent(out)  :: ISHIFTPi(NLINESi)
        character(kind=c_char,len=1),   intent(out) :: IXFIXFPi(NLINESi)
        integer(c_int), intent(out)  :: LINESIZEi(NLINESi)
        character(kind=c_char,len=1),   intent(out) :: AUTOi(NLINESi)
        
        real(c_double), intent(out)  :: OTHER1x(2,NLINESi)
        real(c_double), intent(out)  :: OTHER2x(2,NLINESi)
        real(c_double), intent(out)  :: WLVACi(NLINESi)
        real(c_double), intent(out)  :: CENTERi(NLINESi)
        real(c_double), intent(out)  :: CONCENi(NLINESi)
        integer(c_int), intent(out)  :: NELIONi(NLINESi)
        real(c_float), intent(out)  :: GAMMARi(NLINESi)
        real(c_float), intent(out)  :: GAMMASi(NLINESi)
        real(c_float), intent(out)  :: GAMMAWi(NLINESi)
        character(kind=c_char,len=1),   intent(out) :: REFi(5,NLINESi)
        real(c_float), intent(out)  :: REFx(NLINESi)
        integer(c_int), intent(out) :: NBLOi(NLINESi)
        integer(c_int), intent(out) :: NBUPi(NLINESi)
        integer(c_int), intent(out) :: ISO1i(NLINESi)
        real(c_float), intent(out)  :: X1i(NLINESi)
        integer(c_int), intent(out) :: ISO2i(NLINESi)
        real(c_float), intent(out)  :: X2i(NLINESi)
        real(c_float), intent(out) :: GFLOGi(NLINESi)
        real(c_float), intent(out) :: XJi(NLINESi)
        real(c_float), intent(out) :: XJPi(NLINESi)
        real(c_float), intent(out) :: CODEi(NLINESi)
        real(c_float), intent(out) :: ELOi(NLINESi)
        real(c_float), intent(out) :: GFi(NLINESi)
        real(c_float), intent(out) :: GSi(NLINESi)
        real(c_float), intent(out) :: GRi(NLINESi)
        real(c_float), intent(out) :: GWi(NLINESi)
        real(c_float), intent(out) :: DWLi(NLINESi)
        real(c_float), intent(out) :: DGFLOGi(NLINESi)
        real(c_float), intent(out) :: DGAMMARi(NLINESi)
        real(c_float), intent(out) :: DGAMMASi(NLINESi)
        real(c_float), intent(out) :: DGAMMAWi(NLINESi)
        real(c_float), intent(out) :: DWLISOi(NLINESi)
        integer(c_int), intent(out) :: ISOSHIFTi(NLINESi)
        real(c_float), intent(out) :: EXTRA3i(NLINESi)

        PARAMETER (kw=99)
        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        INTEGER TYPE,J,K
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)
    !   correction 18 May 2011  plus new version of subroutine ionpots.
        COMMON /POTION/POTION(999)
!   C     COMMON /POTION/POTION(594)
        REAL*8 POTION
        DIMENSION CODEX(17)
        DIMENSION DELLIM(7)
        DIMENSION NTENS(10)
        DATA NTENS/1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000/
        DATA CODEX/1.,2.,2.01,6.,6.01,12.,12.01,13.,13.01,14.,14.01,20.,20.01,8.,11.,5.,19./
        DATA DELLIM/100.,30.,10.,3.,1.,.3,.1/

        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        I = 1
        DO WHILE (1.eq.1)
          READ(1,end=100)LINDAT8,LINDAT4
            WLi(I) = WL
            Ei(I) = E
            EPi(I) = EP

            ! Craete label strings
            ! WRITE(SLABEL,'(A8)') LABEL(1)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=1,8
            !   LABELi(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SLABEL,'(A2)') LABEL(2)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=9,10
            !   LABELi(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            LABELx(1,I) = LABEL(1)
            LABELx(2,I) = LABEL(2)


            ! Craete label strings
            ! WRITE(SLABEL,'(A8)') LABELP(1)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=1,8
            !   LABELPi(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SLABEL,'(A2)') LABELP(2)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=9,10
            !   LABELPi(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            LABELPx(1,I) = LABELP(1)
            LABELPx(2,I) = LABELP(2)
            
            READ(COTHER1,'(2I5)')ISHIFT,ISHIFTP
            READ(COTHER2,'(A6,I1,A3)')IXFIXFP,LINESIZE,AUTO

            IF(I.LE.10) WRITE(6,'(I5,I5,A6,I5,A3)') ISHIFT,ISHIFTP,IXFIXFP,LINESIZE,AUTO

            ISHIFTi(I) = ISHIFT
            ISHIFTPi(I) = ISHIFTP
            LINESIZEi(I) = LINESIZE

            ! IXFIXFPi(I) = IXFIXFP(I)
            ! AUTOi(I) = AUTO

            ! WRITE(SLABEL,'(A5)') OTHER1(1)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=1,5
            !   OTHER1i(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SLABEL,'(A5)') OTHER1(2)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=6,10
            !   OTHER1i(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SOTHER1,'(2A5)') OTHER1
            ! SOTHER1 = SOTHER1//c_null_char
            ! DO J=1,10
            !   OTHER1i(J,I) = SOTHER1(J:J)
            ! END DO
            OTHER1x(1,I) = OTHER1(1)
            OTHER1x(2,I) = OTHER1(2)
            
            ! WRITE(SLABEL,'(A7)') OTHER2(1)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=1,7
            !   OTHER2i(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SLABEL,'(A3)') OTHER2(2)
            ! SLABEL = SLABEL//c_null_char
            ! K = 1
            ! DO J=8,10
            !   OTHER2i(J,I) = SLABEL(K:K)
            !   K = K + 1
            ! END DO
            ! WRITE(SOTHER2,'(A6,A3)') OTHER2
            ! SOTHER2 = SOTHER2//c_null_char
            ! DO J=1,10
            !   OTHER2i(J,I) = SOTHER2(J:J)
            ! END DO
            OTHER2x(1,I) = OTHER2(1)
            OTHER2x(2,I) = OTHER2(2)

            WLVACi(I) = WLVAC
            CENTERi(I) = CENTER
            CONCENi(I) = CONCEN
            NELIONi(I) = NELION
            GAMMARi(I) = GAMMAR
            GAMMASi(I) = GAMMAS
            GAMMAWi(I) = GAMMAW

            WRITE(SLABEL,'(A5)') REF

            ! SLABEL = SLABEL//c_null_char
            ! DO J=1,5
            !   REFi(J,I) = SLABEL(J:J)
            ! END DO

            REFx(I) = REF
            NBLOi(I) = NBLO
            NBUPi(I) = NBUP
            ISO1i(I) = ISO1
            X1i(I) = X1
            ISO2i(I) = ISO2
            X2i(I) = X2
            GFLOGi(I) = GFLOG
            XJi(I) = XJ
            XJPi(I) = XJP
            CODEi(I) = CODE
            ELOi(I) = ELO
            GFi(I) = GF
            GSi(I) = GS
            GRi(I) = GR
            GWi(I) = GW
            DWLi(I) = DWL
            DGFLOGi(I) = DGFLOG
            DGAMMARi(I) = DGAMMAR
            DGAMMASi(I) = DGAMMAS
            DGAMMAWi(I) = DGAMMAW
            DWLISOi(I) = DWLISO
            ISOSHIFTi(I) = ISOSHIFT
            EXTRA3i(I) = EXTRA3

            ! IF(I.LE.10) WRITE(6,'(A8,A2)')LABELx(:,I)
            I = I + 1
            ! WRITE(6,*) 'A', OTHER1x(1,I), OTHER1x(2,I)
        END DO
100     CONTINUE
        CLOSE(UNIT=1)
      end subroutine readfile14_20
          
      subroutine writefile14_20(&
        s, NLINESi,&
        WLi,Ei,EPi,LABELi,LABELPi,OTHER1i,OTHER2i,&
        WLVACi,CENTERi,CONCENi,NELIONi,GAMMARi,GAMMASi,GAMMAWi,REFi,&
        NBLOi,NBUPi,ISO1i,X1i,ISO2i,X2i,GFLOGi,XJi,XJPi,CODEi,ELOi,GFi,GSi,GRi,GWi,&
        DWLi,DGFLOGi,DGAMMARi,DGAMMASi,DGAMMAWi,DWLISOi,ISOSHIFTi,EXTRA3i) bind(c, name='writefile14_20')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)
      
        integer(c_int), intent(in), value :: NLINESi
            
        real(c_double), intent(in) :: WLi(NLINESi)
        real(c_double), intent(in) :: Ei(NLINESi)
        real(c_double), intent(in) :: EPi(NLINESi)
        real(c_double), intent(in) :: LABELi(2,NLINESi)
        real(c_double), intent(in) :: LABELPi(2,NLINESi)
        real(c_double), intent(in) :: OTHER1i(2,NLINESi)
        real(c_double), intent(in) :: OTHER2i(2,NLINESi)
        real(c_double), intent(in) :: WLVACi(NLINESi)
        real(c_double), intent(in) :: CENTERi(NLINESi)
        real(c_double), intent(in) :: CONCENi(NLINESi)
        integer(c_int), intent(in) :: NELIONi(NLINESi)
        real(c_float), intent(in) :: GAMMARi(NLINESi)
        real(c_float), intent(in) :: GAMMASi(NLINESi)
        real(c_float), intent(in) :: GAMMAWi(NLINESi)
        real(c_float), intent(in) :: REFi(NLINESi)
        integer(c_int), intent(in) :: NBLOi(NLINESi)
        integer(c_int), intent(in) :: NBUPi(NLINESi)
        integer(c_int), intent(in) :: ISO1i(NLINESi)
        real(c_float), intent(in)  :: X1i(NLINESi)
        integer(c_int), intent(in) :: ISO2i(NLINESi)
        real(c_float), intent(in)  :: X2i(NLINESi)
        real(c_float), intent(in) :: GFLOGi(NLINESi)
        real(c_float), intent(in) :: XJi(NLINESi)
        real(c_float), intent(in) :: XJPi(NLINESi)
        real(c_float), intent(in) :: CODEi(NLINESi)
        real(c_float), intent(in) :: ELOi(NLINESi)
        real(c_float), intent(in) :: GFi(NLINESi)
        real(c_float), intent(in) :: GSi(NLINESi)
        real(c_float), intent(in) :: GRi(NLINESi)
        real(c_float), intent(in) :: GWi(NLINESi)
        real(c_float), intent(in) :: DWLi(NLINESi)
        real(c_float), intent(in) :: DGFLOGi(NLINESi)
        real(c_float), intent(in) :: DGAMMARi(NLINESi)
        real(c_float), intent(in) :: DGAMMASi(NLINESi)
        real(c_float), intent(in) :: DGAMMAWi(NLINESi)
        real(c_float), intent(in) :: DWLISOi(NLINESi)
        integer(c_int), intent(in) :: ISOSHIFTi(NLINESi)
        real(c_float),  intent(in) :: EXTRA3i(NLINESi)

        PARAMETER (kw=99)
        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
               WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
             NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
               DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,DWLISO,ISOSHIFT,EXTRA3
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT4(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT4(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        CHARACTER*10 COTHER1,COTHER2
        EQUIVALENCE (COTHER1,OTHER1(1)),(COTHER2,OTHER2(1))
        INTEGER TYPE,JJ
        EQUIVALENCE (GAMMAS,ASHORE),(GAMMAW,BSHORE)
        EQUIVALENCE (GF,G,CGF),(TYPE,NLAST),(GAMMAR,XSECT,GAUNT)
    !   correction 18 May 2011  plus new version of subroutine ionpots.
        COMMON /POTION/POTION(999)
        REAL*8 POTION
        DIMENSION CODEX(17)
        DIMENSION DELLIM(7)
        DIMENSION NTENS(10)
        DATA NTENS/1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000/
        DATA CODEX/1.,2.,2.01,6.,6.01,12.,12.01,13.,13.01,14.,14.01,20.,20.01,8.,11.,5.,19./
        DATA DELLIM/100.,30.,10.,3.,1.,.3,.1/

        open(UNIT=1,FILE=c_to_f_string(s),STATUS='NEW',FORM='UNFORMATTED')
        JJ = 0
        do I=1,NLINESi

            WL = WLi(I)
            E = Ei(I)
            EP = EPi(I)
            LABEL(1) = LABELi(1,I)
            LABEL(2) = LABELi(2,I)
            LABELP(1) = LABELPi(1,I)
            LABELP(2) = LABELPi(2,I)
            OTHER1(1) = OTHER1i(1,I)
            OTHER1(2) = OTHER1i(2,I)
            OTHER2(1) = OTHER2i(1,I)
            OTHER2(2) = OTHER2i(2,I)

            WLVAC = WLVACi(I)
            CENTER = CENTERi(I)
            CONCEN = CONCENi(I)
            NELION = NELIONi(I)
            GAMMAR = GAMMARi(I)
            GAMMAS = GAMMASi(I)
            GAMMAW = GAMMAWi(I)
            REF = REFi(I)
            NBLO = NBLOi(I)
            NBUP = NBUPi(I)
            ISO1 = ISO1i(I)
            X1 = X1i(I)
            ISO2 = ISO2i(I)
            X2 = X2i(I)
            GFLOG = GFLOGi(I)
            XJ = XJi(I)
            XJP = XJPi(I)
            CODE = CODEi(I)
            ELO = ELOi(I)
            GF = GFi(I)
            GS = GSi(I)
            GR = GRi(I)
            GW = GWi(I)
            DWL = DWLi(I)
            DGFLOG = DGFLOGi(I)
            DGAMMAR = DGAMMARi(I)
            DGAMMAS = DGAMMASi(I)
            DGAMMAW = DGAMMAWi(I)
            DWLISO = DWLISOi(I)
            ISOSHIFT =ISOSHIFTi(I)
            EXTRA3 = EXTRA3i(I)

            ! IF(I.LE.10) WRITE(6,'(A8,A2)') LABEL
            ! WRITE(6,*) 'B',OTHER1
            write(1) LINDAT8,LINDAT4
            JJ = JJ + 1

        end do
        close(unit=1)
    end subroutine writefile14_20

    subroutine readfile93(&
        s, &
        NLINESi,LENGTHi,IFVACi,IFNLTEi,N19i,TURBVi,DECKJi,IFPREDi,&
        WLBEGi,WLENDi,RESOLUi,RATIOi,RATIOLGi,CUTOFFi,LINOUTi) bind(c, name='readfile93')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(out) :: NLINESi
        integer(c_int), intent(out) :: LENGTHi
        integer(c_int), intent(out) :: IFVACi
        integer(c_int), intent(out) :: IFNLTEi
        integer(c_int), intent(out) :: N19i
        real(c_double), intent(out) :: TURBVi
        real(c_double), intent(out) :: DECKJi(7,99)
        integer(c_int), intent(out) :: IFPREDi
        real(c_double), intent(out) :: WLBEGi
        real(c_double), intent(out) :: WLENDi
        real(c_double), intent(out) :: RESOLUi
        real(c_double), intent(out) :: RATIOi
        real(c_double), intent(out) :: RATIOLGi
        real(c_double), intent(out) :: CUTOFFi
        integer(c_int), intent(out) :: LINOUTi

        PARAMETER (kw=99)
        DIMENSION DECKJ(7,kw)
        REAL*8 WLBEG,WLEND,RESOLU,RATIO,RATIOLG
        DATA DECKJ/kw*0.,kw*0.,kw*0.,kw*0.,kw*0.,kw*0.,kw*0./
     
        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        READ(1)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,&
            WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT

        NLINESi  = NLINES
        LENGTHi  = LENGTH
        IFVACi   = IFVAC
        IFNLTEi  = IFNLTE
        N19i     = N19
        TURBVi   = TURBV
        DO J=1,7
            DECKJi(J,:) = DECKJ(J,:)
        END DO
        IFPREDi  = IFPRED
        WLBEGi   = WLBEG
        WLENDi   = WLEND
        RESOLUi  = RESOLU
        RATIOi   = RATIO
        RATIOLGi = RATIOLG
        CUTOFFi  = CUTOFF
        LINOUTi  = LINOUT

    end subroutine readfile93

    subroutine writefile93(&
        s, &
        NLINESi,LENGTHi,IFVACi,IFNLTEi,N19i,TURBVi,DECKJi,IFPREDi,&
        WLBEGi,WLENDi,RESOLUi,RATIOi,RATIOLGi,CUTOFFi,LINOUTi) bind(c, name='writefile93')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        integer(c_int), intent(in) :: NLINESi
        integer(c_int), intent(in) :: LENGTHi
        integer(c_int), intent(in) :: IFVACi
        integer(c_int), intent(in) :: IFNLTEi
        integer(c_int), intent(in) :: N19i
        real(c_double), intent(in) :: TURBVi
        real(c_double), intent(in) :: DECKJi(7,99)
        integer(c_int), intent(in) :: IFPREDi
        real(c_double), intent(in) :: WLBEGi
        real(c_double), intent(in) :: WLENDi
        real(c_double), intent(in) :: RESOLUi
        real(c_double), intent(in) :: RATIOi
        real(c_double), intent(in) :: RATIOLGi
        real(c_double), intent(in) :: CUTOFFi
        integer(c_int), intent(in) :: LINOUTi

        PARAMETER (kw=99)
        DIMENSION DECKJ(7,kw)
        REAL*8 WLBEG,WLEND,RESOLU,RATIO,RATIOLG
        DATA DECKJ/kw*0.,kw*0.,kw*0.,kw*0.,kw*0.,kw*0.,kw*0./
     
        open(UNIT=1,FILE=c_to_f_string(s),STATUS='NEW',FORM='UNFORMATTED')

        NLINES = NLINESi 
        LENGTH = LENGTHi 
        IFVAC  = IFVACi  
        IFNLTE = IFNLTEi 
        N19    = N19i    
        TURBV  = TURBVi  
        DECKJ  = DECKJi
        IFPRED  = IFPREDi  
        WLBEG   = WLBEGi  
        WLEND   = WLENDi  
        RESOLU  = RESOLUi  
        RATIO   = RATIOi   
        RATIOLG = RATIOLGi 
        CUTOFF  = CUTOFFi  
        LINOUT  = LINOUTi  

        WRITE(1)NLINES,LENGTH,IFVAC,IFNLTE,N19,TURBV,DECKJ,IFPRED,&
            WLBEG,WLEND,RESOLU,RATIO,RATIOLG,CUTOFF,LINOUT

        close(unit=1)
    end subroutine writefile93

    subroutine readspecbin_init(s,&
        TEFFi,GLOGi,WBEGINi,RESOLUi,&
        NWLi,IFSURFi,NMUi,NEDGEi,WLEDGEi,&
        NLINESi) bind(c, name='readspecbin_init')    
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)

        real(c_double), intent(out) :: TEFFi
        real(c_double), intent(out) :: GLOGi
        real(c_double), intent(out) :: WBEGINi
        real(c_double), intent(out) :: RESOLUi
        integer(c_int), intent(out) :: NWLi
        real(c_double), intent(out) :: IFSURFi
        real(c_double), intent(out) :: NMUi
        real(c_double), intent(out) :: NEDGEi
        real(c_double), intent(out) :: WLEDGEi
        integer(c_int), intent(out) :: NLINESi

        REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,WLEDGE,RATIO,SWL
        REAL*8 QMU(40),XMU(20),NEDGE,IFSURF,NMU
        REAL*8 NAV,NMU2,FREQTOWAVE
      
        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
                WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
              NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
                DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3&
         ,ALINEC(99)
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN,RESID
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW
        REAL*4 REF,X1,X2,ELO,GF,GS,GR,GW
        REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
        REAL*4 ALINEC
            
        INTEGER IWL, NWL, I, INMU, NLINES, J

        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        read(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,WLEDGE
        ! WRITE(6,*)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,WLEDGE
        TEFFi   = TEFF
        GLOGi   = GLOG
        WBEGINi = WBEGIN
        RESOLUi = RESOLU
        NWLi    = NWL
        IFSURFi = IFSURF
        NMUi    = NMU
        NEDGEi  = NEDGE
        WLEDGEi = WLEDGE

        NAV=IFSURF/10
        IFSURF=IFSURF-NAV*10
        IF(NAV.EQ.0.)NAV=1
        IF(IFSURF.EQ.3)NMU=1
        NMU2=INT(NMU+NMU)
        RATIO=1.D0+1.D0/RESOLU
        
        DO 6 IWL=1,NWL
            READ(1)(QMU(INMU),INMU=1,NMU2)
        6 CONTINUE
        
        READ(1)NLINES
        NLINESi = NLINES
        close(unit=1)

    end subroutine readspecbin_init

    subroutine readspecbin(&
        s, NWLi, NLINESi,&
        wavei, qmu1i, qmu2i, &
        WLi,DWLi,GFLOGi,DGFLOGi,CODEi,&
        GRi,DGAMMARi,GSi,DGAMMASi,GWi,DGAMMAWi,&
        NELIONi,Ei,EPi,NBLOi,NBUPi,ISO1i,ISO2i,X1i,X2i,XJi,XJPi,&
        RESIDi) bind(c, name='readspecbin')
        use iso_c_binding, only: c_double, c_int, c_char, c_null_char, C_LOC, C_PTR, c_float, c_long
        character(kind=c_char,len=1), intent(in) :: s(*)
      
        integer(c_int), intent(in), value :: NWLi
        integer(c_int), intent(in), value :: NLINESi
      
        real(c_double), intent(out) :: wavei(NWLi)
        real(c_double), intent(out) :: qmu1i(NWLi)
        real(c_double), intent(out) :: qmu2i(NWLi)
      
        real(c_double), intent(out) :: WLi(NLINESi)
        real(c_float), intent(out) :: DWLi(NLINESi)
        real(c_float), intent(out) :: GFLOGi(NLINESi)
        real(c_float), intent(out) :: DGFLOGi(NLINESi)
        real(c_float), intent(out) :: CODEi(NLINESi)
        real(c_float), intent(out) :: GRi(NLINESi)
        real(c_float), intent(out) :: DGAMMARi(NLINESi)
        real(c_float), intent(out) :: GSi(NLINESi)
        real(c_float), intent(out) :: DGAMMASi(NLINESi)
        real(c_float), intent(out) :: GWi(NLINESi)
        real(c_float), intent(out) :: DGAMMAWi(NLINESi)

        integer(c_int), intent(out)  :: NELIONi(NLINESi)
        real(c_double), intent(out) :: Ei(NLINESi)
        real(c_double), intent(out) :: EPi(NLINESi)
        integer(c_int), intent(out) :: NBLOi(NLINESi)
        integer(c_int), intent(out) :: NBUPi(NLINESi)
        integer(c_int), intent(out) :: ISO1i(NLINESi)
        integer(c_int), intent(out) :: ISO2i(NLINESi)
        real(c_float), intent(out)  :: X1i(NLINESi)
        real(c_float), intent(out)  :: X2i(NLINESi)
        real(c_float), intent(out) :: XJi(NLINESi)
        real(c_float), intent(out) :: XJPi(NLINESi)

        real(c_double), intent(out) :: RESIDi(NLINESi)
      
        REAL*8 TEFF,GLOG,TITLE(74),WBEGIN,RESOLU,WLEDGE,RATIO,SWL
        REAL*8 QMU(40),XMU(20),NEDGE,IFSURF,NMU
        REAL*8 NAV,NMU2,FREQTOWAVE
      
        COMMON /LINDAT/WL,E,EP,LABEL(2),LABELP(2),OTHER1(2),OTHER2(2),&
                WLVAC,CENTER,CONCEN, NELION,GAMMAR,GAMMAS,GAMMAW,REF,&
              NBLO,NBUP,ISO1,X1,ISO2,X2,GFLOG,XJ,XJP,CODE,ELO,GF,GS,GR,GW,&
                DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3&
         ,ALINEC(99)
        REAL*8 LINDAT8(14)
        REAL*4 LINDAT(28)
        EQUIVALENCE (LINDAT8(1),WL),(LINDAT(1),NELION)
        REAL*8 WL,E,EP,WLVAC,CENTER,CONCEN,RESID
        REAL*8 LABEL,LABELP,OTHER1,OTHER2
        REAL*4 GFLOG,XJ,XJP,CODE,GAMMAR,GAMMAS,GAMMAW
        REAL*4 REF,X1,X2,ELO,GF,GS,GR,GW
        REAL*4 DWL,DGFLOG,DGAMMAR,DGAMMAS,DGAMMAW,EXTRA1,EXTRA2,EXTRA3
        REAL*4 ALINEC
            
        INTEGER IWL, NWL, I, INMU, NLINESO, J
      
        open(UNIT=1,FILE=c_to_f_string(s),STATUS='OLD',FORM='UNFORMATTED',POSITION='REWIND')
        read(1)TEFF,GLOG,TITLE,WBEGIN,RESOLU,NWL,IFSURF,NMU,XMU,NEDGE,WLEDGE
        NAV=IFSURF/10
        IFSURF=IFSURF-NAV*10
        IF(NAV.EQ.0.)NAV=1
        IF(IFSURF.EQ.3)NMU=1
        NMU2=NMU+NMU
        RATIO=1.D0+1.D0/RESOLU
        DO 6 IWL=1,NWL
           SWL=WBEGIN*RATIO**(IWL-1)
           FREQTOWAVE=2.99792458D17/SWL**2
           wavei(IWL) = SWL
           READ(1)(QMU(INMU),INMU=1,NMU2)
        !    IF(IWL.EQ.1)WRITE(6,*) '1',QMU(1)*FREQTOWAVE
        !    IF(IWL.EQ.1)WRITE(6,*) '2',QMU(2)*FREQTOWAVE
           qmu1i(IWL) = QMU(1)*FREQTOWAVE
           qmu2i(IWL) = QMU(2)*FREQTOWAVE
      6 CONTINUE

        READ(1)NLINESO
        
        DO I=1,NLINESO
           READ(1)LINDAT8,LINDAT
      
           WLi(I) = WL
           DWLi(I) = DWL
           GFLOGi(I) = GFLOG
           DGFLOGi(I) = DGFLOG
           CODEi(I) = CODE
           GRi(I) = GR
           DGAMMARi(I) = DGAMMAR
           GSi(I) = GS
           DGAMMASi(I) = DGAMMAS
           GWi(I) = GW
           DGAMMAWi(I) = DGAMMAW

           NELIONi(I) = NELION
           Ei(I) = E
           EPi(I) = EP
           NBLOi(I) = NBLO
           NBUPi(I) = NBUP
           ISO1i(I) = ISO1
           ISO2i(I) = ISO2
           X1i(I) = X1
           X2i(I) = X2
           XJi(I) = XJ
           XJPi(I) = XJP
           
           RESID=CENTER/CONCEN
        !    write(6,*)RESID,CENTER,CONCEN
        !    write(6,*)RESIDi(I)
           RESIDi(I) = RESID
        END DO
      
        CLOSE(UNIT=1)
      
    end subroutine readspecbin

    ! wrap more functions here
    ! ...
    
end module
    