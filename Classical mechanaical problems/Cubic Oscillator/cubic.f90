program cubic
    implicit real*8(a-h,o-z)
    f(x)=-w0**2*x - beta*x**2
    ! Constants
    dt=0.01
    tf=10
    w0=1.0
    beta=0.4
    ! Variable
    t=0.0
    x=1.0
    v=0.0
   
    nsteps=(tf-t)/dt+1

    ! Open file for output data
    open(1,file='cubic.dat')

    do i=1,nsteps
    write(1,*) t,v,x
   
    v=v+f(x)*dt
    x=x+v*dt
    t=t+dt

    enddo
    close(1)
end program cubic
