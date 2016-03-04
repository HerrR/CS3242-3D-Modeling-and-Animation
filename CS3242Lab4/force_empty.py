from math import *

def euler(x, dt, dxdt): 
    x_new = [0, 0, 0]

    # TODO: implement euler method here
    
    return x_new
    
def compute_force_derivatives(force, mass, velocity):
    dxdt = [0, 0, 0]
    dvdt = [0, 0, 0]    

    # TODO: compute location derivative dxdt

    # TODO: compute velocity derivative dvdt

    
    return [dxdt, dvdt]
    
    
def compute_force(const, mass, velocity):
    gravity = const['gravity']
    viscous_drag = const['viscous_drag']
    
    force = [0, 0, 0]

    # TODO: compute total force
    # You can ignore the viscous drag
    
    return force
    
def make_state(x, v):
    return [x, v]
    
def get_location(state):
    return state[0]

def get_velocity(state):
    return state[1]
    
def simulate_force(state, const, mass, dt):
    x = get_location(state)
    v = get_velocity(state)

    # TODO: compute force, derivatives, and make an Euler step on the state
    
    # TODO: make one step Euler
    
    return make_state(x, v)
    
def maya_reset(params):

    initial_location = params['initial_location']
    initial_velocity = params['initial_velocity']
            
    objects = cmds.ls(sl=True)
    if objects == []:
        print('* Please select at least an object.')
        return
        
    for i, o in enumerate(objects):
        
        cmds.select(o)        
        
        if not cmds.attributeQuery('velocityX', n=o, exists=True):
            cmds.addAttr(longName='velocityX', defaultValue=initial_velocity[0])
            cmds.addAttr(longName='velocityY', defaultValue=initial_velocity[1])
            cmds.addAttr(longName='velocityZ', defaultValue=initial_velocity[2])
        else:
            cmds.setAttr(o + '.velocityX', initial_velocity[0])
            cmds.setAttr(o + '.velocityY', initial_velocity[1])
            cmds.setAttr(o + '.velocityZ', initial_velocity[2])
        
        cmds.move(initial_location[0], initial_location[1], initial_location[2])
        
    cmds.select(objects)
    print 'Attributes reset.'
    
def maya_create_destination(const, params):
    # destination is when the object touches the ground again
    # assume no viscous drag
    gravity = const['gravity']
    v = params['initial_velocity']

    # TODO: compute the time the particle touches the ground
    t = 0
    
    x = v[0] * t
    z = v[2] * t
        
    cube = cmds.polyCube(h=2, w=2, d=2)
    cmds.select(cube)
    cmds.move(x, 0.0, z)
    
def maya_move(const, params, time_step):
            
    mass = params['mass']
    initial_velocity = params['initial_velocity']
            
    objects = cmds.ls(sl=True)
    if objects == []:
        print('* Please select at least an object.')
        return
    
    trajectory = cmds.ls('trajectory')
    
    for i, o in enumerate(objects):
        
        cmds.select(o)
        
        x = cmds.getAttr(o + '.translateX')
        y = cmds.getAttr(o + '.translateY')
        z = cmds.getAttr(o + '.translateZ')
        
        if not cmds.attributeQuery('velocityX', n=o, exists=True):
            cmds.addAttr(longName='velocityX', defaultValue=initial_velocity[0])
            cmds.addAttr(longName='velocityY', defaultValue=initial_velocity[1])
            cmds.addAttr(longName='velocityZ', defaultValue=initial_velocity[2])
        
        vx = cmds.getAttr(o + '.velocityX')
        vy = cmds.getAttr(o + '.velocityY')
        vz = cmds.getAttr(o + '.velocityZ')
        
        loc = [x, y, z]
        vel = [vx, vy, vz]
        state = make_state(loc, vel)
        state = simulate_force(state, const, mass, time_step)
        
        old_loc = loc
        loc = get_location(state)
        
        cmds.move(loc[0], loc[1], loc[2])
                
        vel = get_velocity(state)        
        cmds.setAttr(o + '.velocityX', vel[0])
        cmds.setAttr(o + '.velocityY', vel[1])
        cmds.setAttr(o + '.velocityZ', vel[2])
        
        # draw trajectory for the first object
        if i == 0:
            if trajectory == []:
                cv = cmds.curve(point=[old_loc, loc], degree=1)                
                cmds.rename(cv, 'trajectory')
            else:
                cmds.curve('trajectory', point=[loc], degree=1, append=True)
        
    # keep all objects selected
    cmds.select(objects)

# NOTE: remember to call reset every time initial velocity is changed.
const = {'gravity' : 0.98, 'viscous_drag' : 0.0}
params = {'mass' : 1, 'initial_location' : [0, 0, 0], 'initial_velocity' : [5, 5, 0.0]}
maya_create_destination(const, params)
maya_reset(params)

maya_move(const, params, 1.0)
