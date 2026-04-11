program nonlinear_resonance
    implicit real*8(a-h,o-z)

    ! Nonlinear Force: -w0^2*x - 2*gamma*v - beta*x^2 + F0*cos(wd*t)
    f(x,v,t,wd)=-w0**2*x - 2.0d0*gamma*v - beta*x**2 + F0*cos(wd*t)

    ! Constants
    w0 = 1.0d0
    gamma = 0.05d0   ! Low damping to emphasize the tilt
    beta = 0.4d0     ! Nonlinearity coefficient (stiffening spring)
    F0 = 0.2d0       ! Driver strength
    
    dt = 0.02
    tf = 150.0       ! Longer time to ensure steady state
    t_steady = 120.0 

    ! Frequency sweep limits
    w_start = 0.5d0
    w_end = 2.0d0
    dw = 0.01d0
    n_freq = (w_end - w_start)/dw + 1

    ! --- SWEEP UP ---
    open(1,file='sweep_up.dat')
    x = 0.0d0
    v = 0.0d0
    do k=1,n_freq
        wd = w_start + (k-1)*dw
        call get_amplitude(wd, x, v, amp, w0, gamma, beta, F0, dt, tf, t_steady)
        write(1,*) wd, amp
    enddo
    close(1)

    ! --- SWEEP DOWN ---
    open(2,file='sweep_down.dat')
    x = 0.0d0
    v = 0.0d0
    do k=n_freq,1,-1
        wd = w_start + (k-1)*dw
        call get_amplitude(wd, x, v, amp, w0, gamma, beta, F0, dt, tf, t_steady)
        write(2,*) wd, amp
    enddo
    close(2)

end program nonlinear_resonance

! Subroutine to evolve the system and find steady-state amplitude
subroutine get_amplitude(wd, x, v, amp, w0, gamma, beta, F0, dt, tf, t_steady)
    implicit real*8(a-h,o-z)
    t = 0.0d0
    amp = 0.0d0
    nsteps = tf/dt
    
    do i=1,nsteps
        ! Euler-Cromer Integration
        accel = -w0**2*x - 2.0d0*gamma*v - beta*x**2 + F0*cos(wd*t)
        v = v + accel*dt
        x = x + v*dt
        t = t + dt
        
        ! Capture peak amplitude in the steady-state window
        if (t .gt. t_steady) then
            if (dabs(x) .gt. amp) amp = dabs(x)
        endif
    enddo
end subroutine get_amplitude
