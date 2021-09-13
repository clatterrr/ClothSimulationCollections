# ClothSimulationCollections

Large Steps In Cloth Simulation

features:

- compute stretch force,shear force,bend force using triangle element
- conjugate gradient method solve implicit dynamic equation
- only for study,no rendering and optimazation

external library:

freeglut,glew,eigen

substep:

- step01:only shear and stretch force,no damping,explicit step
- step02:shear and stretch force with damping,explicit step
- step03:shear and stretch force with damping,implict step,solve inverse dense matrix
- step04:shear and stretch force with damping,implict step,conjugate gradient to solve sparse matrix
- step05:bend,shear and stretch force with damping,implict step,conjugate gradient to solve sparse matrix

