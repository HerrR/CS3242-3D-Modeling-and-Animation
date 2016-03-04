from maya import cmds
from math import *

def euler(x, dt, dxdt): 
    x_new = [0, 0, 0]

    # TODO: implement the Euler method here

    return x_new

def make_state(x):
    return [x]
    
def get_location(state):
    return state[0]
    
def compute_circle_derivatives(state, angular_velocity):
    x = get_location(state)
    
    dxdt = [0, 0, 0]           

    # TODO: compute the derivative dxdt from the current location and the angular velocity
    
    return dxdt
    
def simulate_circle(state, angular_velocity, dt):
    x = get_location(state)
        
    dxdt = compute_circle_derivatives(state, angular_velocity)
    
    x = euler(x, dt, dxdt)
    
    return make_state(x)
    
def maya_move(angular_velocity, time_step):
            
    objects = cmds.ls(sl=True)
    if objects == []:
        print('* Please select at least an object.')
        return
        
    trajectory = cmds.ls('trajectory')
    
    for i, o in enumerate(objects):
        x = cmds.getAttr(o + '.translateX')
        y = cmds.getAttr(o + '.translateY')
        z = cmds.getAttr(o + '.translateZ')
    
        loc = [x, y, z]
        state = make_state(loc)
                
        state = simulate_circle(state, angular_velocity, time_step)
        
        old_loc = loc
        loc = get_location(state)
        
        cmds.select(o)
        cmds.move(loc[0], loc[1], loc[2])
        
        # draw trajectory for the first object
        if i == 0:
            if trajectory == []:
                cv = cmds.curve(point=[old_loc, loc], degree=1)                
                cmds.rename(cv, 'trajectory')
            else:
                cmds.curve('trajectory', point=[loc], degree=1, append=True)
        
    # keep all objects selected
    cmds.select(objects)
    
maya_move(0.1 * pi, 1.0)
