"Compute inertia matrix for simple solids: cube, sphere and cylinder"
"https://gist.github.com/awesomebytes/39a4ba6c64956a1aa9bd"

def get_cube_inertia_matrix(mass, x, y, z):

    ixx = (1.0 / 12.0) * (y**2 + z**2) * mass
    iyy = (1.0 / 12.0) * (x**2 + z**2) * mass
    izz = (1.0 / 12.0) * (x**2 + y**2) * mass
    ixy = 0.0
    ixz = 0.0
    iyz = 0.0
    return [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]]


def get_sphere_inertia_matrix(mass, radius):
   
    ixx = iyy = izz = (2.0 / 5.0) * radius**2 * mass
    ixy = 0.0
    ixz = 0.0
    iyz = 0.0
    return [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]]


def get_cylinder_inertia_matrix(mass, radius, height):

    ixx = (1.0 / 12.0) * (3.0 * radius**2 + height**2) * mass
    iyy = (1.0 / 12.0) * (3.0 * radius**2 + height**2) * mass
    izz = (1.0 / 2.0) * radius**2 * mass
    ixy = 0.0
    ixz = 0.0
    iyz = 0.0
    return [[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]]

if __name__ == '__main__':
    import sys
    if len(sys.argv) < 2:
        print ("Inertia matrix as:")
        print ("[[ixx, ixy, ixz], [ixy, iyy, iyz], [ixz, iyz, izz]]")
        print ("Cube inertia of cube mass 1.0 and all dimensions 1.0:")
        print (get_cube_inertia_matrix(1.0, 1.0, 1.0, 1.0))
        print ("Sphere inertia of mass 1.0 and radius 1.0:")
        print (get_sphere_inertia_matrix(1.0, 1.0))
        print ("Cylinder inertia of mass 1.0 and radius 1.0 and height 1.0:")
        print (get_cylinder_inertia_matrix(1.0, 1.0, 1.0))
    elif len(sys.argv) == 3:
        print (get_sphere_inertia_matrix(float(sys.argv[1]),
                                  float(sys.argv[2])))
    elif len(sys.argv) == 4:
        print (get_cylinder_inertia_matrix(float(sys.argv[1]),
                                    float(sys.argv[2]),
                                    float(sys.argv[3])))
    elif len(sys.argv) == 5:
        print ("Cube inertia matrix:")
        print (get_cube_inertia_matrix(float(sys.argv[1]),
                                float(sys.argv[2]),
                                float(sys.argv[3]),
                                float(sys.argv[4])))
    else:
        print ("Get the inertia matrix of simple shapes.")
        print ("Usage:")
        print ("For a sphere:")
        print (sys.argv[0] + " mass radius")
        print ("For a cylinder:")
        print (sys.argv[0] + " mass radius height")
        print ("For a cube:")
        print (sys.argv[0] + " mass x y z")
        print ("\nFor example:\n" + sys.argv[0] + "1.0 1.0 1.0 1.0")
        print ("[[0.16666666666666666, 0.0, 0.0], [0.0, 0.16666666666666666, 0.0], [0.0, 0.0, 0.16666666666666666]]")