program tester 

    use mod_global 
    use mod_PGA3D 
    use mod_PGA3Dsupport

    IMPLICIT NONE 

    type(PGA3D_T)           :: mv, mv2, mv3
    integer(kind=irg)       :: i 
    real(kind=dbl)          :: a, b, c, d, x, y, z 

    call PGA3D_initialize(verbose=.TRUE.)

    mv = PGA3D_T()
    call mv%setcomp( -15.D0, 1 )
    call mv%log( pre='test')

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
    mv3 = mv.adds.a 
    call mv3%log( pre='adds')
    write (*,*) 'norm(mv)  = ', mv%norm()
    write (*,*) 'inorm(mv) = ', mv%inorm()
    mv3 = mv%normalized()
    call mv3%log( pre='normalized')
    write (*,*) 'norm(mv3)  = ', mv3%norm()
    mv3 = rotor(cPi/2.0, E12)
    call mv3%log( pre='rotor')
    mv3 = plane (1.D0, 2.D0, 3.D0, 4.D0)
    call mv3%log( pre='plane')
    call getplane(mv3, a, b, c, d)
    write (*,*) 'plane indices : ', a, b, c, d 
    mv3 = point(1.D0, 2.D0, 3.D0)
    call mv3%log( pre='point')
    call getpoint(mv3, x, y, z)
    write (*,*) 'point coordinates : ', x, y, z


end program