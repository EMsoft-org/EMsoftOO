program tester 

    use mod_global 
    use mod_PGA3D 
    use mod_PGA3Dsupport
    use mod_polyhedra
    use mod_STL
    use mod_MCA 

    IMPLICIT NONE 

    type(STL_T)             :: STL 
    type(MCA_T)             :: MCA 
    type(polyhedron_T)      :: shape     
    type(PGA3D_T)           :: mv, mv2, mv3, tr, pt, rot, axz, orig, px, l, p, &
                                rline, rpoint, rplane, pop, pot, poc, cr 
    integer(kind=irg)       :: i, j, k, dims(3), nthr, ntriangles 
    real(kind=dbl)          :: a, b, c, d, x, y, z, dk 
    real(kind=dbl),allocatable            :: gr(:),sf(:,:,:)
    complex(kind=dbl),allocatable         :: shamp(:,:,:)
    real(kind=sgl),allocatable            :: shampreal(:,:,:)
    character(fnlen)        :: sname 
    character(80)           :: header 


    call PGA3D_initialize()
    
    sname = 'snub_cube'
    ! sname = 'cube'
    ! sname = 'icosidodecahedron'
    shape = polyhedron_T( sname, 5.D0)
    ! call shape%polyhedron_info()

    dims = (/ 128, 128, 128 /)
    ! dims = (/ 10, 10, 10 /)
    ! allocate(sf(-dims(1):dims(1),-dims(2):dims(2),-dims(3):dims(3)))
    ! x = 5.D0
    ! call shape%polyhedron_shapefunction(sf, dims, x)

    ! ! do i=-dims(1),dims(1)
    ! !     write (*,"(21F4.1)") sf(i,:,0)
    ! ! end do

    allocate(shamp(-dims(1):dims(1)-1,-dims(2):dims(2)-1,-dims(3):dims(3)-1))
    dk = 0.1D0
    nthr = 8
    call shape%polyhedron_shapeamplitude(shamp, dims, dk, nthr)
    ! do i=-dims(1),dims(1)-1
    !     write (*,*) real(shamp(i,3,5))
    ! end do  

    allocate(shampreal(2*dims(1)+1,2*dims(2)+1,2*dims(3)+1))
    do i=1,2*dims(1)
        do j=1,2*dims(2)
            do k=1,2*dims(3)
                shampreal(i,j,k) = real(shamp(i-dims(1)-1,j-dims(2)-1,k-dims(3)-1))**2
                ! shampreal(i,j,k) = sf(i-dims(1)-1,j-dims(2)-1,k-dims(3)-1)
            end do 
        end do 
    end do 

    MCA = MCA_T()
    call MCA%doMCA( shampreal, 2*dims+1, sngl(dk), 0.01*125.0**2 )

    header = 'Test rendering'
    sname = 'trial.stl' 
    ntriangles = MCA%getNtriangles()
    STL = STL_T(sname, header, ntriangles,MCAlist=MCA%getMCAptr()) 

    ! write (*,*) minval(real(shamp)), maxval(real(shamp))

    ! open(unit=10,file='cube.shamp',status='unknown',form='unformatted')
    ! write (10) real(shamp)
    ! close(unit=10,status='keep')

    ! p = plane(1.D0,1.D0,1.D0,-1.D0)
    ! p = p%normalized()
    ! write (*,*) 'distance to origin = ',p%inorm()


    ! sname = 'cuboctahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'dodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'icosahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'icosidodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'octahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   
    ! sname = 'rhombicosidodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   
    ! sname = 'rhombicuboctahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'rhombitruncated_cuboctahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()    
    ! sname = 'rhombitruncated_icosidodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   
    ! sname = 'snub_cube'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()  
    ! sname = 'snub_dodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   
    ! sname = 'tetrahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()  
    ! sname = 'truncated_cube'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()  
    ! sname = 'truncated_dodecahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   
    ! sname = 'truncated_icosahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()  
    ! sname = 'truncated_octahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()  
    ! sname = 'truncated_tetrahedron'
    ! shape = polyhedron_T( sname )
    ! call shape%polyhedron_info()   



!     mv = PGA3D_T()
!     call mv%setcomp( -15.D0, 1 )
!     call mv%log( pre='test')

!     call mv%log( pre='initial')
!     mv2 = .reverse.mv
!     call mv2%log( pre='reverse')
!     mv2 = .dual.mv 
!     call mv2%log( pre='dual')
!     mv2 = .involute.mv 
!     call mv2%log( pre='involute')
!     mv3 = conjg(mv)
!     call mv3%log( pre='conjg')
!     mv3 = mv * mv2
!     call mv3%log( pre='product')
!     mv3 = mv.wedge.mv2
!     call mv3%log( pre='wedge')
!     mv3 = mv.vee.mv2
!     call mv3%log( pre='vee')
!     mv3 = mv.inner.mv2
!     call mv3%log( pre='inner')
!     a = 20.D0
! ! the following doesn't work yet
!     ! mv3 = a.smul.mv 
!     ! call mv3%log( pre='smul')
!     mv3 = mv.muls.a 
!     call mv3%log( pre='muls')
!     mv3 = mv + mv2
!     call mv3%log( pre='add')
!     mv3 = mv - mv2
!     call mv3%log( pre='subtract')
!     mv3 = mv.adds.a 
!     call mv3%log( pre='adds')
!     write (*,*) 'norm(mv)  = ', mv%norm()
!     write (*,*) 'inorm(mv) = ', mv%inorm()
!     mv3 = mv%normalized()
!     call mv3%log( pre='normalized')
!     write (*,*) 'norm(mv3)  = ', mv3%norm()
    ! rot = rotor(cPi/2.D0, E1*E2)
    ! axz = E1.wedge.E2 
    ! orig = axz.wedge.E3
    ! px = point(0.D0,0.D0,0.D0)
    ! ! line = orig.vee.px 
    ! p = plane(0.D0,0.D0,-1.D0,-3.D0)
    ! mv = p .wedge. px
    ! call mv%log()
    ! mv2 = p.inner.E3
    ! call mv2%log()
    ! write(*,*) mv%getcomp(15)
    ! px = point(0.D0,0.D0, 6.D0)
    ! mv = p .wedge. px
    ! call mv%log()
    ! write(*,*) mv%getcomp(15)

    ! px = point(0.D0,0.D0,0.D0)
    ! ! line = orig.vee.px 
    ! p = plane(0.D0,0.D0,-1.D0, 3.D0)
    ! mv = p .wedge. px
    ! call mv%log()
    ! write(*,*) mv%getcomp(15)
    ! px = point(0.D0,0.D0,-6.D0)
    ! mv = p .wedge. px
    ! call mv%log()
    ! write(*,*) mv%getcomp(15)


    ! rline = rot * line * conjg(rot)
    ! rpoint = rot * px * conjg(rot)
    ! rplane = rot * p * conjg(rot)


    ! call px%log( pre = 'px')
    ! call orig%log( pre = 'orig')
    ! call axz%log( pre = 'axz')
    ! call line%log( pre = 'line')
    ! call p%log( pre = 'plane')
    ! call rot%log( pre = 'rotor')
    ! call rline%log( pre = 'rotated line')
    ! call rpoint%log( pre = 'rotated point')
    ! call rplane%log( pre = 'rotated plane')
    ! pop = (p.inner.px)*p 
    ! call pop%log( pre = 'point on plane (not normalized)')
    ! pop = pop%normalized()
    ! call getpoint(pop, x, y, z)
    ! call pop%log( pre = 'point on plane')
    ! write (*,*) 'point coordinates : ', x, y, z
    ! pot = pointontorus(0.D0, 0.D0)
    ! call pot%log( pre = 'point on torus')
    ! call getpoint(pot, x, y, z)
    ! write (*,*) 'point coordinates : ', x, y, z

    ! tr = translator(5.D0, E1*E0)
    ! do i=0,11
    !     x = dble(i)/12.D0
    !     cr = circle(x, 5.D0, E1*E2)
    !     poc = tr * cr * E123 * conjg(cr) * conjg(tr)
    !     call poc%log()
    !     ! call getpoint(cr, x, y, z)
    !     ! write (*,*) 'point coordinates : ', x, y, z
    ! end do 

    ! p = plane(0.D0, 0.D0, 1.D0, 0.D0)
    ! l = p*E0123
    ! call l%log()

    ! pt = point(0.D0, 0.D0, 5.D0)
    ! p = plane(0.D0, 0.D0, 1.D0, -10.D0)
    ! mv = normalthru(pt, p)
    ! call mv%log()
    ! mv = projectonto(pt, p)
    ! call mv%log()

    ! pt = point(0.D0, 0.D0, 0.D0) 
    ! tr = point(0.D0, 0.D0, 5.D0) 
    ! write (*,*) 'distance = ', distpoints(pt, tr)
    ! l = pt.vee.tr 
    ! write(*,*) 'norm of line ', l%norm()

    ! mv = plane(0.D0, 1.D0, 1.D0, 5.D0)
    ! mv2= plane(0.D0, 1.D0, 0.D0, -20.D0)
    ! write (*,*) 'distance = ', distplanes(mv, mv2)
    ! write (*,*) 'angle     = ', angleplanes(mv, mv2)/dtor
    ! write (*,*) 'angle     = ', angleplaneline(mv, l)/dtor
    ! write (*,*) 'distance  = ', distplaneline(mv2, l)
    ! write (*,*) 'oriented distance = ', ordisttoplane(tr, mv2)
    ! pt = point(2.D0, 2.D0, 0.D0) 
    ! write (*,*) 'oriented distance = ', ordisttoline(pt, l)

    ! pt = point(0.D0, 0.D0, 0.D0)
    ! tr = translator(5.D0, E1*E0)
    ! pt = tr * pt * conjg(tr)
    ! tr = translator(-2.D0, E2*E0) 
    ! pt = tr * pt * conjg(tr)
    ! tr = translator(4.D0, E3*E0)
    ! mv = tr * pt * conjg(tr)
    ! call getpoint(mv, x, y, z)
    ! write (*,*) 'translated point coordinates : ', x, y, z


    ! mv3 = rotor(cPi/2.0, E12)
    ! call mv3%log( pre='rotor')
    ! mv3 = plane (1.D0, 2.D0, 3.D0, 4.D0)
    ! call mv3%log( pre='plane')
    ! call getplane(mv3, a, b, c, d)
    ! write (*,*) 'plane indices : ', a, b, c, d 
    ! mv3 = point(1.D0, 2.D0, 3.D0)
    ! call mv3%log( pre='point')
    ! call getpoint(mv3, x, y, z)
    ! write (*,*) 'point coordinates : ', x, y, z


end program