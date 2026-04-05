program coupled_oscillator 
    implicit real*8(a-h,o-z)
    f(t,x,y,vx,vy)= -w**2*y-2*w**2*x
    g(t,x,y,vx,vy)= -w**2*x-2*w**2*y
    !parameters
    x=0.5
    y=-0.5
    vx=0.0  
    vy=0.0
    w=1.0
    dt=0.01
    tmax=20.0   
    nsteps=int(tmax/dt)
    t=0.0
    
    open(1,file='coupled_oscillator.dat')
    
    do n=1,nsteps 
        !RK-4 method
        ak1x=vx*dt
        ak1y=vy*dt
        ak1vx=f(t,x,y,vx,vy)*dt
        ak1vy=g(t,x,y,vx,vy)*dt

        ak2x=(vx+0.5*ak1vx)*dt
        ak2y=(vy+0.5*ak1vy)*dt
        ak2vx=f(t+0.5*dt,x+0.5*ak1x,y+0.5*ak1y,vx+0.5*ak1vx,vy+0.5*ak1vy)*dt
        ak2vy=g(t+0.5*dt,x+0.5*ak1x,y+0.5*ak1y,vx+0.5*ak1vx,vy+0.5*ak1vy)*dt

        ak3x=(vx+0.5*ak2vx)*dt
        ak3y=(vy+0.5*ak2vy)*dt
        ak3vx=f(t+0.5*dt,x+0.5*ak2x,y+0.5*ak2y,vx+0.5*ak2vx,vy+0.5*ak2vy)*dt
        ak3vy=g(t+0.5*dt,x+0.5*ak2x,y+0.5*ak2y,vx+0.5*ak2vx,vy+0.5*ak2vy)*dt 

        ak4x=(vx+ak3vx)*dt
        ak4y=(vy+ak3vy)*dt
        ak4vx=f(t+dt,x+ak3x,y+ak3y,vx+ak3vx,vy+ak3vy)*dt
        ak4vy=g(t+dt,x+ak3x,y+ak3y,vx+ak3vx,vy+ak3vy)*dt

        x=x+(ak1x+2.0*ak2x+2.0*ak3x+ak4x)/6.0
        y=y+(ak1y+2.0*ak2y+2.0*ak3y+ak4y)/6.0
        vx=vx+(ak1vx+2.0*ak2vx+2.0*ak3vx+ak4vx)/6.0
        vy=vy+(ak1vy+2.0*ak2vy+2.0*ak3vy+ak4vy)/6.0
        t=t+dt
        
        write(1,100) t,x,y,vx,vy
    enddo
    
    close(1)
    
    100 format(5E20.10)
    
end program coupled_oscillator

