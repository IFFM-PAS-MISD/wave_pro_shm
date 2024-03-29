Source files and their organization
***********************************

Source files are located in ``src`` folder and corresponding subfolders:

* ``wave_pro_shm.mlapp`` file contains source related to gui prepared in MATLAB app designer.
* ``mesher`` folder contains utilities for meshing.
* ``wave_pro_shm`` folder contains subfolders:

  * ``solid_3d_engine`` folder contains main functions behind the time domain spectral element method implementation.
 


.. automodule:: src

.. autoclass:: wave_pro_shm

.. automodule:: src.wave_pro_shm


Solid 3D Engine
===============

List of functions in the ``solid3d_engine`` folder.

.. automodule:: src.wave_pro_shm.solid3d_engine

.. autofunction:: det_jacp
.. autofunction:: dofs3d
.. autofunction:: dofs3dfifi
.. autofunction:: elasticity_matrix_isotropy
.. autofunction:: elasticity_matrix_piezo
.. autofunction:: excitation_signals_at_actuators
.. autofunction:: hanning_signal
.. autofunction:: inv_jacp
.. autofunction:: piezoelectricity_matrices
.. autofunction:: pzt_coupling_matrices
.. autofunction:: save_pvd_animation
.. autofunction:: save_vtu_frame
.. autofunction:: shape3d_prim_coeff
.. autofunction:: shape3d_prim_coeff_single_node
.. autofunction:: solid3d_engine
.. autofunction:: sparse_block_diagonal_shape_derivatives
.. autofunction:: strains_spec_p
.. autofunction:: stresses_spec_p

Utilities for meshing
=====================

List of functions in the ``mesher`` folder.

.. automodule:: src.mesher

.. autofunction:: element2d_numbering
.. autofunction:: element3d_gll_nodes_vector_method
.. autofunction:: element3d_numbering
.. autofunction:: gll
.. autofunction:: gmsh2sem_mesh
.. autofunction:: local_global_nodes_map
.. autofunction:: local_global_nodes_map2dofs
.. autofunction:: locate_elements_in_physical_entities
.. autofunction:: nodal_base_change_3d
.. autofunction:: read_msh_mesh
.. autofunction:: retrieve_mesh_damage_elements
.. autofunction:: retrieve_mesh_transducer_elements
.. autofunction:: sem2paraview_mesh
.. autofunction:: split_interface_nodes
.. autofunction:: vandermonde
