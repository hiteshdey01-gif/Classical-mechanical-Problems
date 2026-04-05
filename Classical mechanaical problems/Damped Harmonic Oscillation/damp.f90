
    implicit real*8(a-h,o-z)
    f(x,v)=-w0**2*x-g*v
    open(1, file='overdamp.dat', status='unknown')
    w0=1.0
    g=3.0
    !write(*,*)'INPUT ffinal time'
    !read(*,*)tf
    tf=20
    dt= 0.001
    nstep=int(tf/dt)+1
    t=0.0
    x=1.0
    v=0.0
    do i=1,nstep
        write(1,*)t,x,v
        v=v+f(x,v)*dt
        x=x+v*dt
        t=t+dt
    enddo
    stop 

end
