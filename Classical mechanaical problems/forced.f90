program forced
implicit real*8(a-h,o-z)
    f(x,v,t)= -(k/m)*x-g*v+f0*cos(w*t)
    open(1,file="forced.dat")

    m = 1.0          !Mass
    k = 1.0
    g = 0.5          ! damping coefficient 
    f0 = 7.0         ! driving force amplitude
    w= 3.0           ! driving frequency

    dt = 0.001       ! Time Interval
    tf = 40
    nstep = int(tf/dt) + 1


    x = 1.0
    v = 0.0
    t = 0.0

    do i = 1, nstep
        v = v + f(x,v,t)*dt
        x = x + v*dt
        t = t + dt

    write(1,*) t, x, v
      
    enddo
    stop
end program forced