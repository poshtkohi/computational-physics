!    character(*) meta ...
!    print *, meta

!	All rights reserved to Alireza Poshtkohi (c) 2017-2023.
!	Email: arp@poshtkohi.info
!	Website: http://www.poshtkohi.info
! ----------------------------------------
program main
    common /com2/p(2000000,5),ip(2000000)
    common /elec/cn(200),cnd(200)
    common /field/qhf(200)
    common /mesh/ns,nd,nx,nx1
    common /dope/dn1
    common /bias/vd
    common /size/dx,xs,xd,xl
    common /time/dt,t
    common /jtime/jt,jtl
    common /imp1/cimp
    common /num/inum
    common /mpi/my_rank,num_procs
    !common /time/tmax
    !real*8 dt,t

    integer ierr, my_rank, num_procs, real_procs, d_one_process, global_inum
    integer max_particles_per_process, max_excess_inum_particles_per_process, remaining_particles
    integer i, j, k, process_index, process_inum, particle_index
    !real particles
    real, allocatable::p_master(:)
    real, allocatable::p_slave(:)
    dimension cn_receive(200)

    include 'mpif.h'

    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD, my_rank, ierr)
    call mpi_comm_size(MPI_COMM_WORLD, num_procs, ierr)

    !if((my_rank.eq.0).and.(num_procs.lt.2)) then
    !    print *, "Number of MPI process must be at least 2"
    !    call mpi_finalize(ierr)
    !    stop
    !endif

    max_excess_inum_particles_per_process = 1;
    my_inum = 0
    d_one_process = 0
    i = 0
    j = 0
    k = 0

    ns = 30
    nd = 80
    nx = 110
    dn1 = 2.e+23
    dx = 50.e-10
    dt = 2.e-15
    tem = 300
    vd = 0.5
    cimp = 1.e22
    tmax = 2.e-9
    jtl = 500

    if(num_procs.eq.1) then
        ! Initialization
        call config
        call param
        call initia
        !call log_particles("serial-initia.txt")
        ! EMC simulation
        call cpu_time(start_t)

        do jt=1,jtl
            t = dt * float(jt)
            call emcd
            call renew
            !print *, "inum", inum
            !stop
            call charge
            call poisso
            !print *,"t",t,"jt",jt
        end do

        call cpu_time(end_t)
        write(*, *) "Elapsed CPU time in seconds = ", end_t - start_t

        call output

        call mpi_finalize(ierr)
        stop
        return
    endif

    !write(*, *) "header size ", sizeof(h)
    !write(*, *) "t ", sizeof(t)

    ! print *, dn1

    ! Initialization
    call config
    call param

    print *, "main at rank ", my_rank

    if(my_rank.eq.0) then
        !-------- Initialization --------!
        call initia
        ! First distribute the particle across slaves
        write(*, *) "inum", inum
        real_procs = num_procs - 1
        max_particles_per_process = inum/real_procs
        remaining_particles = mod(inum, real_procs)
        ! Broadcast max_particles_per_process
        call mpi_bcast(max_particles_per_process, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        ! Allocate partciles_master and p_slave. First sizeof(real) is for inum
        d_one_process = sizeof(real) + (max_particles_per_process + max_excess_inum_particles_per_process) * 6 * sizeof(real)
        allocate(p_master(d_one_process * num_procs), stat=ierr)
        allocate(p_slave(d_one_process), stat=ierr)
        ! Send particle information to slaves
        ! First assign the particles in a round-robin fashion to slaves
        ! The indices of i and j are used for particle and processes respectively.
        j = 2
        do i=1,inum
            process_index = 1 + (6 * (max_particles_per_process + max_excess_inum_particles_per_process) + 1) * (j - 1) ! ignores the first process
            process_inum = ifix(p_master(process_index))
            p_master(process_index + 6 * process_inum + 1) = p(i, 1)
            p_master(process_index + 6 * process_inum + 2) = p(i, 2)
            p_master(process_index + 6 * process_inum + 3 ) = p(i, 3)
            p_master(process_index + 6 * process_inum + 4) = p(i, 4)
            p_master(process_index + 6 * process_inum + 5) = p(i, 5)
            p_master(process_index + 6 * process_inum + 6) = ip(i)
            p_master(process_index) = p_master(process_index) + 1
            !write(*, *) "at rank", my_rank, "process_index", process_index, "process_inum", p_master(process_index), "j", j
            if(num_procs.ge.3) then
                if(j.eq.num_procs) then
                    j = 2
                else
                    j = j + 1
                endif
            endif
        end do
        ! Print the information to check their soundness
        !call log_particles('serial.txt')
        ! Scatter the information to slaves
        call mpi_scatter(p_master, d_one_process, MPI_BYTE, p_slave, d_one_process, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
        ! Prepares data allocated arrays for the EMC simulation
        deallocate(p_master, stat=ierr)
        deallocate(p_slave, stat=ierr)
        d_one_process = sizeof(real) + (max_particles_per_process + max_excess_inum_particles_per_process) * 2 * sizeof(real)
        allocate(p_master(d_one_process * num_procs), stat=ierr)
        allocate(p_slave(d_one_process), stat=ierr)
        !-------- EMC Simulation --------!
        ! EMC simulation
        call cpu_time(start_t)
        ! Main loop
        do 20 jt=1,jtl
            ! Generate and send a new t to all salves
            t = dt * float(jt)
            call mpi_bcast(t, 1, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
            k = 0
            !call MPI_Reduce(k, global_inum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            !write(*, *) "at rank", my_rank, "global_inum", global_inum
            !go to 20
            ! Receives particle information after ecmd/renew on slaves
            call mpi_gather(p_slave, d_one_process, MPI_BYTE, p_master, d_one_process, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
            inum = 0
            k = 1
            do i=2,num_procs
                process_index = 1 + (2 * (max_particles_per_process + max_excess_inum_particles_per_process) + 1) * (i - 1) ! ignores the first process
                my_inum = ifix(p_master(process_index))
                do j=1,my_inum
                    !write(*, *) "at rank", my_rank, "process_index", process_index, "my_inum", my_inum, "k", k
                    p(k, 5) = p_master(2 * (j - 1) + 1 + process_index)
                    ip(k) = ifix(p_master(2 * (j - 1) + 2 + process_index))
                    k = k + 1
                end do
                inum = inum + my_inum
                !write(*, *) "at rank", my_rank, "my_inum", my_inum
            end do
            if(inum.eq.0) exit
            !write(*, *) "at rank", my_rank, "inum", inum
            ! Print the information to check their soundness
            !call log_particles('serial.txt')
            !go to 30
            ! Executes charge and poisso equations
            !call charge
            ! First update charges
            do i=1,200
                cn(i) = 0.
            end do
            call mpi_reduce(cn, cn_receive, 200, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            do i=1,200
                cn(i) = cn_receive(i)
            end do
            call poisso
            call mpi_bcast(qhf, 200, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
        20 continue
        ! Compute elapsed time
        call cpu_time(end_t)
        write(*, *) "Elapsed CPU time in seconds = ", end_t - start_t
        ! Terminate the slaves
        if(inum.eq.0) go to 100
        t = 0.
        call mpi_bcast(t, 1, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
        ! Output to file
        call output
        ! Free allocated arrays
        deallocate(p_master, stat=ierr)
        deallocate(p_slave, stat=ierr)
    else
        !-------- Initialization --------!
        ! Receives max_particles_per_process
        call mpi_bcast(max_particles_per_process, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
        !go to 30
        ! Receives local particle information from master
        d_one_process = sizeof(real) + (max_particles_per_process + max_excess_inum_particles_per_process) * 6 * sizeof(real)
        allocate(p_slave(d_one_process), stat=ierr)
        call mpi_scatter(0, 0, 0, p_slave, d_one_process, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
        ! Copy the information to local particles
        process_inum = ifix(p_slave(1))
        inum = process_inum
        !write(*, *) "at rank", my_rank, "p_slave(1)", p_slave(1), "max_particles_per_process", max_particles_per_process
        do i=1,process_inum
            particle_index = 6 * (i - 1) + 1 ! ignores the first index
            p(i, 1) = p_slave(particle_index + 1)
            p(i, 2) = p_slave(particle_index + 2)
            p(i, 3) = p_slave(particle_index + 3)
            p(i, 4) = p_slave(particle_index + 4)
            p(i, 5) = p_slave(particle_index + 5)
            ip(i) = ifix(p_slave(particle_index + 6))
        end do
        ! Print the information to check their soundness
        !call log_particles('parallel.txt')
        ! Prepares data allocated arrays for the EMC simulation
        deallocate(p_slave, stat=ierr)
        d_one_process = sizeof(real) + (max_particles_per_process + max_excess_inum_particles_per_process) * 2 * sizeof(real)
        allocate(p_slave(d_one_process), stat=ierr)
        !-------- EMC Simulation --------!
        do 40
            ! Receives a new t
            call mpi_bcast(t, 1, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
            if(t.eq.0.) exit
            ! Execute ecmd and renew routines on each slave
            !write(*, *) "inum before", inum
            call emcd
            call renew
            call charge
            !call MPI_Reduce(inum, global_inum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            !go to 40
            !write(*, *) "inum after", inum
            ! Send back updated particle information to the master
            p_slave(1) = inum
            do i=1,inum
                p_slave(2 * (i - 1) + 1 + 1) = p(i, 5) ! 1 is for inum index
                p_slave(2 * (i - 1) + 2 + 1) = float(ip(i))
            end do
            call mpi_gather(p_slave, d_one_process, MPI_BYTE, 0, 0, 0, 0, MPI_COMM_WORLD, ierr)
            ! Reduce cn all slaves to master
            call mpi_reduce(cn, cn_receive, 200, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
            call mpi_bcast(qhf, 200, MPI_FLOAT, 0, MPI_COMM_WORLD, ierr)
            if(inum.eq.0) go to 100
            ! Print the information to check their soundness
            !call log_particles('parallel.txt')
            !go to 30
        40 continue
        ! Free allocated arrays
        deallocate(p_slave, stat=ierr)
    endif

    !30 continue
    !write(*, *) "at rank", my_rank, "t", t, "max_particles_per_process", max_particles_per_process

    call mpi_barrier(MPI_COMM_WORLD, ierr)

    !write(*, *) "t",
    100 continue
    !call MPI_Abort(MPI_COMM_WORLD, 0)
    call mpi_finalize(ierr)
    !call exit(-1)
    stop

end program main
! ----------------------------------------
! Logs particle information to a file
subroutine log_particles(filename)
    common /com2/p(2000000,5),ip(2000000)
    common /num/inum
    character(*) :: filename

    open(10, file=filename)

    do i=1,inum
        do j=1,5
            write(10, *) "i", i, "j", j , "val", p(i, j)
        end do
        write(10, *) "i", i, "j", 6 , "val", ip(i)
    end do

    flush(10)
    close(10)

    return

end subroutine log_particles
! ----------------------------------------
! Configuration of device
subroutine config
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    common /ptcl/np1,epp ! np1 is number of super particles.
    common /dope/dn1
    common /elec/cn(200),cnd(200)
    common /mtrix1/a1(200),b1(200),c1(200)
    common /mtrix2/bt(200),gm(200)
    common /mpi/my_rank,num_procs

    parameter(q = 1.60219e-19)

    print *, "config at rank ", my_rank

    nx1 = nx + 1
    np1 = 2000!50000 !100, 2000
    epp = (dn1 * dx)/float(np1)
    xs = dx * (float(ns) - 0.5)
    xd = dx * (float(nd) - 1.5)
    xl = dx * float(nx)

    ! Doping at the cathode & anode
    do i=1,nx1
        if((i.le.ns).or.(i.ge.nd)) then
            cnd(i) = dn1
        end if
    end do

    ! Coefficients for Poisson equation
    do i=1,nx1
        a1(i) = 1.
        b1(i) = -2
        c1(i) = 1.
    end do

    ! LU decomposition of a1, b1 and c1
    bt(2) = b1(2)
    gm(2) = c1(2) / bt(2)
    do i=3,nx
        bt(i) = b1(i) - a1(i) * gm(i - 1)
        gm(i) = c1(i) / bt(i)
    end do

    return

end subroutine config
! ----------------------------------------
! Physical constants & parameters
subroutine param
    common /rate/hwo,hwij,hwe,de,swk(2,8,1000)
    common /gamma/gm
    common /ek/smh(2),hhm(2),ec(2)
    common /nonp/af(2),af2(2),af4(2)
    common /hm2/hm(2)
    common /dri/qh
    common /temp/tem
    common /bf/bktq
    common /imp1/cimp
    common /imp2/dq2
    common /emax/iemax
    common /random/iseed
    common /mpi/my_rank,num_procs
    real no,nij,ne

    parameter(pi = 3.14159, q=1.60219e-19)
    parameter(h = 1.05459e-34, bk = 1.38066e-23)
    parameter(ep0 = 8.85419e-12, am0= 9.10953e-31)

    print *, "param at rank ", my_rank

    ! initial value for random number generator
    iseed = 38476

    ! Effetive masses
    am1 = 0.067 * am0 ! gamma
    am2 = 0.350 * am0 ! L

    ! Energy minima of bands
    ec(1) = 0.
    ec(2) = 0.29

    ! Energy gap & dieletric constant
    eg = 1.424
    eps = 12.90 * ep0
    epf = 10.92 * ep0
    ep = 1./(1./epf-1./eps)

    ! Parameters for phonon scatterings
    rou = 5360.
    sv = 5240.
    cl = rou * sv * sv
    z2 = 4.
    z3 = 3.
    da = 7. * q
    dij = 1.e11 * q
    deq = 1.e11 * q
    hwo = 0.03536
    hwij = 0.03
    hwe = hwij

    ! Non-parabolicity of Gamma & L & X bands
    af(1) = (1.-am1/am0)**2/eg
    af2(1) = 2. * af(1)
    af4(1) = 4. * af(1)
    af(2) = (1.-am2/am0)**2/(eg+ec(2))
    af2(2) = 2. * af(2)
    af4(2) = 4. * af(2)

    bktq = bk*tem/q
    qh = q/h
    smh(1) = sqrt(2.*am1)* sqrt(q)/h
    smh(2) = sqrt(2.*am2)*sqrt(q)/h
    hhm(1) = h/am1/q*h/2.
    hhm(2) = h/am2/q*h/2.
    hm(1) = h/am1
    hm(2)= h/am2

    wo = hwo*q/h
    wij = hwij*q/h
    we = hwe*q/h

    no = 1./(exp(hwo/bktq)-1.)
    nij = 1./(exp(hwij/bktq)-1.)
    ne = 1./(exp(hwe/bktq)-1.)

    dos1 = (sqrt(2.*am1)*sqrt(q)/h)**3/4./pi/pi
    dos2 = (sqrt(2.*am2)*sqrt(q)/h)**3/4./pi/pi

    poe = q/8./pi/ep*q*wo*(no+1.)
    poa = poe*no/(1.+no)
    aco = 2.*pi*da/q*da*bktq/h*q/cl
    ope = pi*dij/wij*dij/rou/q*(nij+1.)
    opa = ope*nij/(1.+nij)
    eqe = pi*deq/we*deq/rou/q*(ne+1.)
    eqa = eqe*ne/(1.+ne)

    ! Parameters for impurity scattering
    qd = sqrt(q*cimp/bktq/eps)
    qd2 = qd*qd
    bimp = 2.*pi*cimp*q*q/h*q/eps/eps

    ! Calculation of scattering rates
    de = 0.002
    iemax = 1000

    do ie=1,iemax
        ei = de*float(ie)
        sei = sqrt(ei)
        ! --------- Gamma-valleys ----------------
        ! ---- Polar optical phonon ----
        ef = ei-hwo
        if(ef.gt.0.) then
            sef=sqrt(ef)
            qmax=sef+sei
            qmin=sei-sef
            swk(1,1,ie)=poe*smh(1)*sei/ei/q*alog(qmax/qmin)
        else
            swk(1,1,ie) = 0.
        endif

        ef = ei+hwo
        sef = sqrt(ef)
        qmax = sef+sei
        qmin = sef-sei
        swk(1,2,ie) = swk(1,1,ie)+poa*smh(1)*sei/ei/q*alog(qmax/qmin)
        ! ---- Non-polar optical phonon ----
        ef = ei-hwij+ec(1)-ec(2)
        if(ef.gt.0.) then
            sef=sqrt(ef*(1.+af(2)*ef))
            swk(1,3,ie)=swk(1,2,ie)+z2*ope*sef*dos2*(1.+2.*af(2)*ef)
        else
            swk(1,3,ie)=swk(1,2,ie)
        endif

        ef=ei+hwij+ec(1)-ec(2)
        if(ef.gt.0.) then
            sef=sqrt(ef*(1.+af(2)*ef))
            swk(1,4,ie)=swk(1,3,ie)+z2*opa*sef*dos2*(1.+2.*af(2)*ef)
        else
            swk(1,4,ie)=swk(1,3,ie)
        endif
        ! ---- Acoustic phonon ----
        ef=ei
        sef=sqrt(ef*(1.+af(1)*ef))
        swk(1,5,ie)=swk(1,4,ie)+aco*sef*dos1*(1.+2.*af(1)*ef)
        ! ---- Impurity scattering ----
        ef=ei
        sef=sqrt(ef*(1.+af(1)*ef))
        ak=smh(1)*sef
        qq=qd2*(4.*ak*ak+qd2)
        wk=bimp/qq*sef*dos1*(1.+2.*af(1)*ef)
        if(wk.gt.1.e14) wk=1.e14
        swk(1,6,ie)=swk(1,5,ie)+wk
        !print *, "ak", ak, "qd2", qd2, "qq", qq, "wk", wk
        ! --------- L-valleys ----------------
        ! ---- Polar optical phonon ----
        ef=ei-hwo
        if(ef.gt.0.) then
            sef=sqrt(ef)
            qmax=sef+sei
            qmin=sei-sef
            swk(2,1,ie)=poe*smh(2)*sei/ei/q*alog(qmax/qmin)
        else
            swk(2,1,ie) = 0.
        endif

        ef=ei+hwo
        sef=sqrt(ef)
        qmax=sef+sei
        qmin=sef-sei
        swk(2,2,ie)=swk(2,1,ie)+poa*smh(2)*sei/ei/q*alog(qmax/qmin)
        ! ---- Non-polar optical phonon ----
        ef=ei-hwe
        if(ef.gt.0.) then
            sef = sqrt(ef*(1.+af(2)*ef))
            swk(2,3,ie)=swk(2,2,ie)+(z2-1.)*eqe*sef*dos2*(1.+2.*af(2)*ef)
        else
            swk(2,3,ie)=swk(2,2,ie)
        endif

        ef=ei+hwe
        if(ef.gt.0.) then
            sef = sqrt(ef*(1.+af(2)*ef))
            swk(2,4,ie)=swk(2,3,ie)+(z2-1.)*eqa*sef*dos2*(1.+2.*af(2)*ef)
        else
            swk(2,4,ie)=swk(2,3,ie)
        endif

        ef=ei-hwij+ec(2)-ec(1)
        if(ef.gt.0.) then
            sef = sqrt(ef*(1.+af(1)*ef))
            swk(2,5,ie)=swk(2,4,ie)+ope*sef*dos1*(1.+2.*af(1)*ef)
        else
            swk(2,5,ie)=swk(2,4,ie)
        endif

        ef=ei+hwij+ec(2)-ec(1)
        if(ef.gt.0.) then
            sef = sqrt(ef*(1.+af(1)*ef))
            swk(2,6,ie)=swk(2,5,ie)+opa*sef*dos1*(1.+2.*af(1)*ef)
        else
            swk(2,6,ie)=swk(2,5,ie)
        endif

        ! ---- Acoustic phonon ----
        ef=ei
        sef=sqrt(ef*(1.+af(2)*ef))
        swk(2,7,ie)=swk(2,6,ie)+aco*sef*dos2*(1.+2.*af(2)*ef)
        ! ---- Impurity scattering ----
        ef=ei
        sef=sqrt(ef*(1.+af(2)*ef))
        ak=smh(2)*sef
        qq=qd2*(4.*ak*ak+qd2)
        wk=bimp/qq*sef*dos2*(1.+2.*af(2)*ef)
        if(wk.gt.1.e14) wk=1.e14
        swk(2,8,ie)=swk(2,7,ie)+wk
    end do

    !print *,"swk(1,2,1)=",swk(1,2,2)
    !stop

    ! --------- Evaluation of Gamma ----------------
    gm=swk(1,6,1)

    do ie=1,iemax
        if(swk(1,6,ie).gt.gm) gm=swk(1,6,ie)
        if(swk(2,8,ie).gt.gm) gm=swk(2,8,ie)
    end do

    do i=1,6
        do ie=1,iemax
            swk(1,i,ie)=swk(1,i,ie)/gm
            !print *,"swk(1,",i,",",ie,")=",swk(1,i,ie)
        end do
    end do

    do i=1,8
        do ie=1,iemax
            swk(2,i,ie)=swk(2,i,ie)/gm
            !print *,"swk(2,",i,",",ie,")=",swk(2,i,ie)
        end do
    end do

    return

end subroutine param
! ----------------------------------------
! Initial condition for particles
subroutine initia
    common /com2/p(2000000,5),ip(2000000)
    common /ek/smh(2),hhm(2),ec(2)
    common /nonp/af(2),af2(2),af4(2)
    common /bf/bktq
    common /gamma/gm
    common /num/inum
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    common /ptcl/np1,epp
    common /elec/cn(200),cnd(200)
    common /field/qhf(200)
    common /mpi/my_rank,num_procs

    parameter(pi = 3.1415927, q=1.60219e-19)

    print *, "initia at rank", my_rank

    n = 0

    do i=1,nx1
        npi = ifix(cnd(i)*dx/epp+0.5)
        if((i.eq.1).or.(i.eq.nx1)) npi = npi/2
        if(npi.eq.0) exit
        do m=1,npi
            n = n + 1
            if(n.gt.2000000) then
                print *,"Number of particles exceeds 2000000", "n",n
                stop
            end if
            iv = 1
            ei = -bktq*alog(rnd())*1.5
            ak = smh(iv)*sqrt(ei*(1.+af(iv)*ei))
            cb = 1.-2.*rnd()
            sb = sqrt(1.-cb*cb)
            fai = 2.*pi*rnd()
            sf = sin(fai)
            cf = cos(fai)
            p(n,1) = ak*cb*sf
            p(n,2) = ak*sb*sf
            p(n,3) = ak*cf
            p(n,4) = -alog(rnd())/gm
            p(n,5) = dx*(rnd()+float(i)-1.5)
            if(i.eq.1) p(n,5) = dx*0.5*rnd()
            if(i.eq.nx1) p(n,5) = xl-dx*0.5*rnd()
            ip(n) = iv
            !print *, "n", n, "ip", ip(n), "p(n,1)", p(n,1)
        end do
        !print *, "i", i
    end do

    inum = n
    do n=inum+1,2000000
        ip(n) = 9
    end do

    do i=1,nx1
        qhf(i) = 0.
    end do

    print *, "n", n, "inum", inum, "epp", epp

    return

end subroutine initia
! ----------------------------------------
! Particle motion during dt
subroutine emcd
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /com2/p(2000000,5),ip(2000000)
    common /gamma/gm
    common /num/inum
    common /time/dt,t
    real kx,ky,kz

    !print *, "ecmd"

    !print *, "inum", inum

    tdt = t + dt
    do 30 n=1,inum
        kx = p(n,1)
        ky = p(n,2)
        kz = p(n,3)
        ts = p(n,4)
        x = p(n,5)
        iv = ip(n)

        t1 = t
     10 if(ts.gt.tdt) go to 20
        tau = ts - t1
        call drift(tau)
        call scat
        t1 = ts
        ts = t1 - alog(rnd())/gm
        go to 10
     20 tau = tdt - t1
        call drift(tau)

        p(n,1) = kx
        p(n,2) = ky
        p(n,3) = kz
        p(n,4) = ts
        p(n,5) = x
        ip(n) = iv

        !print *, "n", n
    30 continue

    return

end subroutine emcd
! ----------------------------------------
! Calculation of drift process
subroutine drift(tau)
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /field/qhf(200)
    common /nonp/af(2),af2(2),af4(2)
    common /ek/smh(2),hhm(2),ec(2)
    common /hm2/hm(2)
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    real kx,ky,kz

    !print *, "drift"

    if(iv.eq.9) return

    i = max(1,min(ifix(x/dx+1.5),nx1))
    dkx = qhf(i)*tau
    hmt = hm(iv)*tau
    sk = kx*kx+ky*ky+kz*kz
    gk = hhm(iv)*sk
    sq = sqrt(1.+af4(iv)*gk)
    x = x + hmt*(kx+0.5*dkx)/sq
    kx = kx + dkx

    call surf

    return

end subroutine drift
! ----------------------------------------
! Boundary condition for particle
subroutine surf
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /size/dx,xs,xd,xl
    real kx,ky,kz

    !print *, "surf"

    if((x.le.0.).or.(x.ge.xl)) then
        iv = 9
    end if

    return

end subroutine surf
! ----------------------------------------
! Calculation of scattering process
subroutine scat
    common /rate/hwo,hwij,hwe,de,swk(2,8,1000)
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /ek/smh(2),hhm(2),ec(2)
    common /nonp/af(2),af2(2),af4(2)
    common /imp2/qd2
    common /emax/iemax
    common /size/dx,xs,xd,xl
    real kx,ky,kz,ki,kf

    parameter(pi = 3.14159)

    !print *, "scat"

    ei=e
    if(ei.eq.0.) return

    ki=sqrt(sk)
    ie=ifix(ei/de)+1
    if(ie.gt.iemax) ie=iemax

    if(iv.eq.1) then
        goto 1000
    elseif(iv.eq.2) then
        goto 2000
    elseif(iv.eq.9) then
        return
    endif
    return

    ! Selection of scattering process
    1000 r1=rnd()

    if(r1.le.swk(1,1,ie)) then
        ef=ei-hwo
        if(ef.le.0.) return
        goto 20
    elseif(r1.le.swk(1,2,ie)) then
        ef=ei+hwo
        goto 20
    elseif(r1.le.swk(1,3,ie)) then
        ef=ei-hwij+ec(1)-ec(2)
        if(ef.le.0.) return
        iv=2
        goto 40
    elseif(r1.le.swk(1,4,ie)) then
        ef=ei+hwij+ec(1)-ec(2)
        if(ef.le.0.) return
        iv=2
        goto 40
    elseif(r1.le.swk(1,5,ie)) then
        ef=ei
        kf=ki
        goto 40
    elseif(r1.le.swk(1,6,ie)) then
        if((x.lt.xs).or.(x.gt.xd)) then
            ef=ei
            r2=rnd()
            cb=1.-r2/(0.5+(1.-r2)*sk/qd2)
            kf=ki
            goto 30
        endif
    endif
    return

    2000 r1=rnd()

    if(r1.le.swk(2,1,ie)) then
        ef=ei-hwo
        if(ef.le.0.) return
        goto 20
    elseif(r1.le.swk(2,2,ie)) then
        ef=ei+hwo
        goto 20
    elseif(r1.le.swk(2,3,ie)) then
        ef=ei-hwe
        if(ef.le.0.) return
        goto 40
    elseif(r1.le.swk(2,4,ie)) then
        ef=ei+hwe
        goto 40
    elseif(r1.le.swk(2,5,ie)) then
        ef=ei-hwij+ec(2)-ec(1)
        if(ef.le.0.) return
        iv=1
        goto 40
    elseif(r1.le.swk(2,6,ie)) then
        ef=ei+hwij+ec(2)-ec(1)
        if(ef.le.0.) return
        iv=1
        goto 40
    elseif(r1.le.swk(2,7,ie)) then
        ef=ei
        kf=ki
        goto 40
    elseif(r1.le.swk(2,8,ie)) then
        ef=ei
        r2=rnd()
        cb=1.-r2/(0.5+(1.-r2)*sk/qd2)
        kf=ki
        goto 30
    endif
    return

    ! Determination of final states

    20 kf=smh(iv)*sqrt(ef*(1.+af(iv)*ef))
    f=2.*ki*kf/(ki-kf)/(ki-kf)
    if(f.le.0.) return
    cb=(1.+f-(1.+2.*f)**rnd())/f

    30 sb=sqrt(1.-cb*cb)
    fai=2.*pi*rnd()
    cf=cos(fai)
    sf=sin(fai)
    skk=sqrt(kx*kx+ky*ky)
    a11=ky/skk
    a12=kx*kz/skk/ki
    a13=kx/ki
    a21=-kx/skk
    a22=ky*kz/skk/ki
    a23=ky/ki
    a32=-skk/ki
    a33=kz/ki
    x1=kf*sb*cf
    x2=kf*sb*sf
    x3=kf*cb
    kx=a11*x1+a12*x2+a13*x3
    ky=a21*x1+a22*x2+a23*x3
    kz=a32*x2+a33*x3
    e=ef
    return

    40 kf=smh(iv)*sqrt(ef*(1.+af(iv)*ef))
    cs=1.-2*rnd()
    sn=sqrt(1.-cs*cs)
    fai=2.*pi*rnd()
    kx=kf*cs
    ky=kf*sn*cos(fai)
    kz=kf*sn*sin(fai)
    e=ef

    return

end subroutine scat
! ----------------------------------------
! Renewing process
subroutine renew
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /com2/p(2000000,5),ip(2000000)
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    common /ptcl/np1,epp
    common /num/inum
    common /mpi/my_rank,num_procs

    !integer npt
    dimension npt(2)
    real kx,ky,kz
    !integer ic

    !print *, "renew"

    np12 = np1/2

    npt(1) = 0
    npt(2) = 0
    ic = 0

    !print *, "inum", inum
    nn = inum
    do 50 n=1,nn
        10 if(ip(inum).eq.9) then
            inum = inum - 1
            go to 10
        endif
        if(n.ge.inum) then
            go to 30
        endif

        if(ip(n).eq.9) then
            do 20 j=1,5
                p(n,j) = p(inum,j)
            20 continue
            ip(n) = ip(inum)
            ip(inum) = 9
            inum = inum - 1
        endif

        30 if(ip(inum).eq.9) then
            inum = inum - 1
            go to 30
        endif
        if(n.ge.inum) then
            go to 60
        endif

        i = max(1,min(ifix(p(n,5)/dx+1.5),nx1))
        if(i.eq.1) then
            ic = 1
        elseif(i.eq.nx1) then
            ic = 2
        else
            go to 50
        endif
        if(npt(ic).ge.np12) then
            do 40 j=1,5
                p(n,j) = p(inum,j)
            40 continue
            ip(n) = ip(inum)
            ip(inum) = 9
            inum = inum - 1
        else
            npt(ic) = npt(ic) + 1
        endif
    50 continue

    60 if(num_procs.gt.1) then
        return
    endif

    !print *, "npt(1)", npt(1), "npt(2)", npt(2), "np12", np12
    do 80 ic=1,2
        ni = np12 - npt(ic)
        !print *, "ni", ni
        if(ni.lt.0) go to 80
        do 70 j=1,ni
            n = inum + j
            !!print *, "n", n
            call create(ic)
            p(n,1) = kx;
            p(n,2) = ky
            p(n,3) = kz
            p(n,4) = ts
            p(n,5) = x
            ip(n) = iv
        70 continue
        inum = inum + ni
    80 continue

    return

end subroutine renew
! ----------------------------------------
! Creating a particle
subroutine create(ic)
    common /com1/kx,ky,kz,sk,ts,x,iv,e
    common /bf/bktq
    common /ek/smh(2),hhm(2),ec(2)
    common /nonp/af(2),af2(2),af4(2)
    common /gamma/gm
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    common /time/dt,t
    real kx,ky,kz

    parameter(pi = 3.14159)

    !print *, "create"

    iv = 1
    ei = -bktq*alog(rnd())*1.5
    ak = smh(iv)*sqrt(ei*(1.+af(iv)*ei))
    cb = rnd()
    sb = sqrt(1.-cb*cb)
    fai = 2.*pi*rnd()
    sf = sin(fai)
    cf = cos(fai)
    kx = ak*cb
    ky = ak*sb*sf
    kz = ak*sb*cf
    ts = t - alog(rnd())/gm

    if(ic.eq.1) then
        x = dx*0.5*rnd()
    elseif(ic.eq.2) then
        kx = -kx
        x = xl - dx*0.5*rnd()
    endif

    return

end subroutine create
! ----------------------------------------
! Charge distribution
subroutine charge
    common /com2/p(2000000,5),ip(2000000)
    common /elec/cn(200),cnd(200)
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    common /num/inum
    common /ptcl/np1,epp
    common /jtime/jt,jtl

    !print *, "charge"

    do 10 i=1,nx1
        cn(i) = 0.
    10 continue

    ! Cloud-in-cell scheme
    do 20 n=1,inum
        if(ip(n).eq.9) go to 20
        x = p(n,5)/dx
        i = max(1,min(ifix(x)+1,nx))
        xb = float(i-1)
        x2 = 1.-(x-xb)
        cn(i) = cn(i) + x2
        cn(i+1) = cn(i + 1) + (1.-x2)
    20 continue

    do 30 i=1,nx1
        cn(i) = cn(i)*epp/dx
        if((i.eq.1).or.(i.eq.nx1)) cn(i) = cn(i)*2.
    30 continue

    return

end subroutine charge
! ----------------------------------------
! Potential calculation
subroutine poisso
    common /elec/cn(200),cnd(200)
    common /field/qhf(200)
    common /pot1/psi(200)
    common /mtrix1/a1(200),b1(200),c1(200)
    common /mtrix2/bt(200),gm(200)
    common /pot2/qeps
    common /dri/qh
    common /bias/vd
    common /mesh/ns,nd,nx,nx1
    common /size/dx,xs,xd,xl
    dimension z(200),h(200)

    !print *, "poisso"

    psi(1) = 0.
    psi(nx1) = vd

    do 10 i=2,nx
        h(i) = (cn(i)-cnd(i))*qeps
    10 continue
    h(2) = h(2)-a1(2)*psi(1)
    h(nx) = h(nx) - c1(nx)*psi(nx1)

    z(2) = h(2)/b1(2)
    do 20 i=3,nx
        z(i) = (h(i)-a1(i)*z(i-1))/bt(i)
    20 continue

    psi(nx) = z(nx)
    do 30 i=nx-1,2,-1
        psi(i) = z(i)-gm(i)*psi(i+1)
    30 continue

    do 40 i=2,nx
        qhf(i) = qh*(psi(i+1)-psi(i-1))/dx/2
    40 continue
    qhf(1) = qhf(2)
    qhf(nx1) = qhf(nx)

    return

end subroutine poisso
! ----------------------------------------
! Output routine
subroutine output
    common /field/qhf(200)
    common /pot1/psi(200)
    common /size/dx,xs,xd,xl
    common /mesh/ns,nd,nx,nx1

    ! Open files
    open(8, file = 'potential.txt')
    open(9, file = 'field.txt')

    do i=2,nx
        write(8,*) i*dx," ",psi(i)
        write(9,*) i*dx," ",qhf(i)
    end do

    flush(8)
    flush(9)

    close(8)
    close(10)

    return
end
! ----------------------------------------
! Random number generator
function rnd()
    common /random/iseed
    parameter(mi = 1048576, in = 1027)
    iseed = mod(in*iseed, mi)
    !rnd = float(iseed)/float(mi)
    rnd = rand(iseed)
    return
end function rnd
! ----------------------------------------

