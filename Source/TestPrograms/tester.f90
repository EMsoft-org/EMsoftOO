program tester 

    use mod_global 
    use mod_PGA3D 

    IMPLICIT NONE 

    type(PGA3D_T)           :: mv, mv2, mv3 
    integer(kind=irg)       :: i 
    real(kind=dbl)          :: a 

    mv = PGA3D_T()

    do i=0,15,2
      call mv%setcomp( dble(i), i)
    end do 
    call mv%setcomp( -15.D0, 4 )
    do i=0,15 
        write (*,*) mv%getcomp(i)
    end do 

    call mv%log( pre='initial')
    mv2 = .reverse.mv
    call mv2%log( pre='reverse')
    mv2 = .dual.mv 
    call mv2%log( pre='dual')
    mv2 = .involute.mv 
    call mv2%log( pre='involute')
    mv3 = conjg(mv)
    call mv3%log( pre='conjg')
    mv3 = mv * mv2
    call mv3%log( pre='product')
    mv3 = mv.wedge.mv2
    call mv3%log( pre='wedge')
    mv3 = mv.vee.mv2
    call mv3%log( pre='vee')
    mv3 = mv.inner.mv2
    call mv3%log( pre='inner')
    a = 20.D0
! the following doesn't work yet
    ! mv3 = a.smul.mv 
    ! call mv3%log( pre='smul')
    mv3 = mv.muls.a 
    call mv3%log( pre='muls')
    mv3 = mv + mv2
    call mv3%log( pre='add')
    mv3 = mv - mv2
    call mv3%log( pre='subtract')
    mv3 = mv + a 
    call mv3%log( pre='adds')



end program