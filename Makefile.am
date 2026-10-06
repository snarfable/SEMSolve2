LDFLAGS = -llapack -lxml2
bin_PROGRAMS = semsolve mk_rand_nodes mk_rand_sphere mk_sphere_hull mk_mesh smooth_mesh smooth_dist convert_mesh testing
semsolve_SOURCES = semsolve.cpp
semsolve_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
mk_rand_nodes_SOURCES = mk_rand_nodes.cpp
mk_rand_nodes_LDADD = libsemsolve.a
mk_rand_sphere_SOURCES = mk_rand_sphere.cpp
mk_rand_sphere_LDADD = libsemsolve.a
mk_sphere_hull_SOURCES = mk_sphere_hull.cpp
mk_sphere_hull_LDADD = libsemlinalg.a
smooth_mesh_SOURCES = smooth_mesh.cpp
smooth_mesh_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
smooth_dist_SOURCES = smooth_dist.cpp
smooth_dist_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
mk_mesh_SOURCES = mk_mesh.cpp
mk_mesh_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
convert_mesh_SOURCES = convert_mesh.cpp
convert_mesh_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
testing_SOURCES = testing.cpp
testing_LDADD = libsemsolve.a libsemlinalg.a -llapack -lblas
lib_LIBRARIES = libsemsolve.a libsemlinalg.a
libsemsolve_a_SOURCES =  basis.cpp basis.h\
			 mesh.cpp mesh.h\
			 ltgmap.cpp ltgmap.h\
			 combinadic.cpp combinadic.h\
			 semerror.cpp semerror.h\
			 semconfig.cpp semconfig.h\
			 integrals.cpp integrals.h\
			 potential.cpp potential.h\
			 geometry.cpp geometry.h\
			 parsing.y lexing.l
libsemlinalg_a_SOURCES = dlinalg.h dlinalg_util.cpp dlinalg_eigen.cpp dlinalg_lnsys.cpp\
			 slinalg.h slinalg_util.cpp slinalg_eigen.cpp slinalg_lnsys.cpp\
			 dlinalg_asm.s
