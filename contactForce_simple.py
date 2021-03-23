'''
This function is from moorepants GaitAnalysisToolkit. Slightly adapted with predefined physical constants.
'''

from sympy import symbols, exp

def contact_force(point, ground, origin):
    """Returns a contact force vector acting on the given point made of
    friction along the contact surface and elastic force in the vertical
    direction.

    Parameters
    ==========
    point : sympy.physics.mechanics.Point
        The point which the contact force should be computed for.
    ground : sympy.physics.mechanics.ReferenceFrame
        A reference frame which represents the inerital ground in 3D space.
        The x axis defines the ground line and positive z is up.
    origin : sympy.physics.mechanics.Point
        An origin point located on the ground line.

    Returns
    =======
    force : sympy.physics.mechanics.Vector
        The contact force between the point and the ground.

    """
    # This is the "height" of the point above the ground, where a negative
    # value means that the point is below the ground.
    z_location = point.pos_from(origin).dot(ground.z)
    

    # The penetration into the ground is mathematically defined as:
    #
    #               { 0 if z_location > 0
    # deformation = {
    #               { abs(z_location) if z_location < 0
    #

    penetration = (abs(z_location) - z_location) / 2

    velocity = point.vel(ground)

    # The addition of "- y_location" here adds a small linear term to the
    # cubic stiffness and creates a light attractive force torwards the
    # ground. This is in place to ensure that gradients can be computed for
    # the optimization used in Ackermann and van den Bogert 2010.

    contact_stiffness = 5.0e7
    contact_damping = 0.85
    contact_friction_coefficient = 1
    friction_scaling_factor = 1
#    contact_friction_coefficient, friction_scaling_factor = \
#        symbols('mu, vs')
 
    vertical_force = (contact_stiffness * penetration ** 3 - z_location) * \
        (1 - contact_damping * velocity.dot(ground.z))
    
    friction_y = -vertical_force * ((2 / (1 + exp(-velocity.dot(ground.y)))) - 1)

    friction_x = -vertical_force * ((2 / (1 + exp(-velocity.dot(ground.x)))) - 1)

    return friction_x * ground.x + friction_y *  ground.y + vertical_force * ground.z
    #return vertical_force * ground.z
