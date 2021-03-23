from sympy import symbols, Eq, solve, pprint, simplify, trigsimp
from sympy.physics.mechanics import ReferenceFrame, dynamicsymbols, Point, inertia, RigidBody, Particle, KanesMethod, mprint
from numpy import linspace, deg2rad, pi, ones, array, sin, set_printoptions
from pydy.system import System
from pydy.viz import Plane, Cylinder, Sphere, VisualizationFrame, Scene
import matplotlib.pyplot as plt
from contactForce_simple import contact_force




N = ReferenceFrame('N')
B = ReferenceFrame('B') 
URA = ReferenceFrame('URA')
LRA = ReferenceFrame('LRA')

ULA = ReferenceFrame('ULA')
LLA = ReferenceFrame('LLA')

URL = ReferenceFrame('URL')
LRL = ReferenceFrame('LRL')

ULL = ReferenceFrame('ULL')
LLL = ReferenceFrame('LLL')

q1, q2, q3 = dynamicsymbols('q1:4')
theta1 = dynamicsymbols('theta1')
omega1 = dynamicsymbols('omega1')


B.orient(N, 'Body', (q1, q2, q3), 'xyz')

URA.orient(B, 'Body', (pi/2-deg2rad(+20), +pi/2, 0), 'yzy')
LRA.orient(URA, 'Axis', (-pi/2+theta1, URA.x))


ULA.orient(B, 'Body', (pi/2-deg2rad(-20), +pi/2, 0), 'yzy')
LLA.orient(ULA, 'Axis', (-pi/2, ULA.x))

URL.orient(B, 'Body', (pi/2-deg2rad(+20), +pi/2, 0), 'yzy')
LRL.orient(URL, 'Axis', (+pi/2, URL.x))

ULL.orient(B, 'Body', (pi/2-deg2rad(-20), +pi/2, 0), 'yzy')
LLL.orient(ULL, 'Axis', (+pi/2, ULL.x))

O = Point('O') 
O.set_vel(N, 0)


G = O.locatenew('G', -6.5 * N.z)

x, y, z = dynamicsymbols('x, y, z')
C = O.locatenew('C', x * N.x + y * N.y + z * N.z)

RA = C.locatenew('RA', 6.0 * B.x + -4 * B.y)
LA = C.locatenew('LA', 6.0 * B.x + 4 * B.y)
RL = C.locatenew('RL', -6.0 * B.x + -4 * B.y)
LL = C.locatenew('LL', -6.0 * B.x + 4 * B.y)

C_URA = RA.locatenew('C_URA', -2.5 * URA.y)
C_ULA = LA.locatenew('C_ULA', -2.5 * ULA.y)
C_URL = RL.locatenew('C_URL', -2.5 * URL.y)
C_ULL = LL.locatenew('C_ULL', -2.5 * ULL.y)

elbow_RA = RA.locatenew('elbow_RA', -5 * URA.y)
elbow_LA = LA.locatenew('elbow_LA', -5 * ULA.y)
knee_RL = RL.locatenew('knee_RL', -5 * URL.y)
knee_LL = LL.locatenew('knee_LL', -5 * ULL.y)

C_LRA = elbow_RA.locatenew('C_LRA', -2.5 * LRA.y)
C_LLA = elbow_LA.locatenew('C_LLA', -2.5 * LLA.y)
C_LRL = knee_RL.locatenew('C_LRL', -2.5 * LRL.y)
C_LLL = knee_LL.locatenew('C_LLL', -2.5 * LLL.y)

paw_RA = elbow_RA.locatenew('paw_RA', -5 * LRA.y)
paw_LA = elbow_LA.locatenew('paw_LA', -5 * LLA.y)
paw_RL = knee_RL.locatenew('paw_RL', -5 * LRL.y)
paw_LL = knee_LL.locatenew('paw_LL', -5 * LLL.y)

ux = dynamicsymbols('u_x')
uy = dynamicsymbols('u_y')
uz = dynamicsymbols('u_z')
u1, u2, u3 = dynamicsymbols('u_1:4')

z1 = Eq(ux, x.diff())
z2 = Eq(uy, y.diff())
z3 = Eq(uz, z.diff())
z4 = Eq(u1, q1.diff())
z5 = Eq(u2, q2.diff())
z6 = Eq(u3, q3.diff())
z7 = Eq(omega1, theta1.diff())

u = solve([z1, z2, z3, z4, z5, z6, z7], x.diff(), y.diff(), z.diff(), q1.diff(), q2.diff(), q3.diff())
mprint(u)

C.set_vel(N, C.pos_from(O).dt(N).subs(u))

B.set_ang_vel(N, B.ang_vel_in(N).subs(u))
URA.set_ang_vel(B, 0)
ULA.set_ang_vel(B, 0)
URL.set_ang_vel(B, 0)
ULL.set_ang_vel(B, 0)

LRA.set_ang_vel(URA, omega1 * URA.x)
LLA.set_ang_vel(ULA, 0)
LRL.set_ang_vel(URL, 0)
LLL.set_ang_vel(ULL, 0)

C_URA.v2pt_theory(RA, N, URA)
C_ULA.v2pt_theory(LA, N, ULA)
C_URL.v2pt_theory(RL, N, URL)
C_ULL.v2pt_theory(LL, N, ULL)

elbow_RA.v2pt_theory(RA, N, URA)
elbow_LA.v2pt_theory(LA, N, ULA)
knee_RL.v2pt_theory(RL, N, URL)
knee_LL.v2pt_theory(LL, N, ULL)

C_LRA.v2pt_theory(elbow_RA, N, LRA)
C_LLA.v2pt_theory(elbow_LA, N, LLA)
C_LRL.v2pt_theory(knee_RL, N, LRL)
C_LLL.v2pt_theory(knee_LL, N, LLL)

paw_RA.v2pt_theory(elbow_RA, N, LRA)
paw_LA.v2pt_theory(elbow_LA, N, LLA)
paw_RL.v2pt_theory(knee_RL, N, LRL)
paw_LL.v2pt_theory(knee_LL, N, LLL)



m, m_link = symbols('m, m_link') # Nybble mass

Ix, Iy, Iz = symbols('I_x, I_y, I_z') # principal moments of inertia

I = inertia(B, Ix, Iy, Iz) # inertia dyadic


Fz_mag = dynamicsymbols('Fmag_z')
g = symbols('g')

Fz = Fz_mag * N.z * g

kdes = [z1.rhs - z1.lhs,
        z2.rhs - z2.lhs,
        z3.rhs - z3.lhs,
        z4.rhs - z4.lhs,
        z5.rhs - z5.lhs,
        z6.rhs - z6.lhs,
        z7.rhs - z7.lhs,
        ]

bodies = []
bodies.append(RigidBody('body', C, B, m, (I, C)))
bodies.append(RigidBody('upper_arm_r', C_URA, URA, m_link, (I,C_URA)))
bodies.append(RigidBody('lower_arm_r', C_LRA, LRA, m_link, (I,C_LRA)))
bodies.append(RigidBody('upper_arm_l', C_ULA, ULA, m_link, (I,C_URA)))
bodies.append(RigidBody('lower_arm_l', C_LLA, LLA, m_link, (I,C_LRA)))
bodies.append(RigidBody('upper_leg_r', C_URL, URL, m_link, (I,C_URL)))
bodies.append(RigidBody('lower_leg_r', C_LRL, LRL, m_link, (I,C_LRL)))
bodies.append(RigidBody('upper_leg_l', C_ULL, ULL, m_link, (I,C_ULL)))
bodies.append(RigidBody('lower_leg_l', C_LLL, LLL, m_link, (I,C_LLL)))


loads = [
        (C, Fz),
        (C_URA, Fz),
        (C_LRA, Fz),
        (C_ULA, Fz),
        (C_LLA, Fz),
        (C_URL, Fz),
        (C_LRL, Fz),
        (C_ULL, Fz),
        (C_LLL, Fz),
        (paw_RA, contact_force(paw_RA, N, G)),
        (paw_LA, contact_force(paw_LA, N, G)),
        (paw_LL, contact_force(paw_LL, N, G)),
        (paw_RL, contact_force(paw_RL, N, G)),
        (LRA, -1000 * omega1 * URA.x),
        ]

kane = KanesMethod(N, (x, y, z, q1, q2, q3, theta1), (ux, uy, uz, u1, u2, u3, omega1), kd_eqs=kdes)
fr, frstar = kane.kanes_equations(bodies, loads=loads)


sys = System(
        kane,
        constants = {
            Ix: 0.1083,
            Iy: 0.1083,
            Iz: 0.1083,
            m: 7,
            m_link: 1,
            g: -9.81,
        },

        times = linspace(0.0, 3, num=90),



        specifieds = {
            Fz_mag : 1.0,
            #theta1: deg2rad(-20)
        })

sys.generate_ode_function(generator='cython') # Speed up integration with Cython

states = []

sys.initial_conditions = {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            q1: 0.0,
            q2: 0.0,
            q3: 0.0,
            ux: 0.0,
            uy: 0.0,
            uz: 0.0,
            u1: 0.0,
            u2: 0.0,
            u3: 0.0,
            theta1: deg2rad(-10),
            omega1: 0.0,
        }

states.append(sys.integrate())
import numpy
numpy.set_printoptions(threshold=numpy.inf)
print(states)

#for i in range(0,3):
#    sys.initial_conditions = {
#                x: states[0+i],
#                y: states[1+i],
#                z: states[2+i],
#                q1: states[3+i],
#                q2: states[4+i],
#                q3: states[5+i],
#                ux: states[6+i],
#                uy: states[7+i],
#                uz: states[8+i],
#                u1: states[9+i],
#                u2: states[10+i],
#                u3: states[11+i],
#            }


#    states.append(sys.integrate())
#    mprint(states[i])

fig, ax = plt.subplots()
ax.plot(sys.times, states[0])
#ax.plot(states)
ax.set_xlabel('time t [s]', fontsize=8)
ax.set_ylabel('pos. q [m] / vel. u_xyz [ms-1] / angular vel. u_123 [s-1]', fontsize=8)
ax.legend(['$q_1$', '$q_2$', '$q_3$', '$u_x$', '$u_y$', '$u_z$', '$u_1$', '$u_2$', '$u_3$'], fontsize=8)
plt.show()


body_shape = Plane(12, 8, color='blue')
link_shape = Cylinder(radius=0.08, length= 5,  color='black')
joint_shape = Sphere(color='black', radius=0.2)
ground_shape = Plane(60, 60, color='white')

viz_objects = []
viz_objects.append(VisualizationFrame('Body_m', B, C, body_shape))

viz_objects.append(VisualizationFrame('right_shouler', B, RA, joint_shape))
viz_objects.append(VisualizationFrame('upper_right_arm', URA, C_URA , link_shape))
viz_objects.append(VisualizationFrame('elbow_right_arm', URA, elbow_RA, joint_shape))
viz_objects.append(VisualizationFrame('lower_right_arm', LRA, C_LRA, link_shape))
viz_objects.append(VisualizationFrame('paw_right_arm', LRA, paw_RA, joint_shape))

viz_objects.append(VisualizationFrame('left_shouler', B, LA, joint_shape))
viz_objects.append(VisualizationFrame('upper_left_arm', ULA, C_ULA , link_shape))
viz_objects.append(VisualizationFrame('elbow_left_arm', ULA, elbow_LA, joint_shape))
viz_objects.append(VisualizationFrame('lower_left_arm', LLA, C_LLA, link_shape))
viz_objects.append(VisualizationFrame('paw_left_arm', LLA, paw_LA, joint_shape))

viz_objects.append(VisualizationFrame('right_hip', B, RL, joint_shape))
viz_objects.append(VisualizationFrame('upper_right_leg', URL, C_URL, link_shape))
viz_objects.append(VisualizationFrame('knee_rigt_leg', URL, knee_RL, joint_shape))
viz_objects.append(VisualizationFrame('lower_right_leg', LRL, C_LRL, link_shape))
viz_objects.append(VisualizationFrame('paw_right_leg', LRL, paw_RL, joint_shape))

viz_objects.append(VisualizationFrame('left_hip', B, LL, joint_shape))
viz_objects.append(VisualizationFrame('upper_left_leg', ULL, C_ULL, link_shape))
viz_objects.append(VisualizationFrame('knee_left_leg', ULL, knee_LL, joint_shape))
viz_objects.append(VisualizationFrame('lower_left_leg', LLL, C_LLL, link_shape))
viz_objects.append(VisualizationFrame('paw_left_leg', LLL, paw_LL, joint_shape))

viz_objects.append(VisualizationFrame('ground', N, G, ground_shape))


scene = Scene(N, O, system=sys)
scene.visualization_frames = viz_objects
scene.display()
