					/* Project node m in local plane g onto each axis in g */
// 					for(p=1;p<i;p++)
// 					{
// 						vec1[p-1] = 0;
// 						for(q=0;q<dim;q++)
// 							vec1[p-1] += ( nodes[ map_in->npmap[g][m] * dim + q ] - stdvec[cmb[0]*dim+q] )
// 									* ( stdvec[cmb[p]*dim+q] - stdvec[cmb[0]*dim+q] );
// 						sum = 0.0;
// 						for(q=0;q<dim;q++)
// 							sum += powf( stdvec[cmb[p]*dim+q] - stdvec[cmb[0]*dim+q], 2.0 );
// 						vec1[p-1] /= sqrt( sum );
// 					}

						/* Project node n in local plane h onto each axis in h */
// 						for(p=1;p<i;p++)
// 						{
// 							vec2[p-1] = 0;
// 							for(q=0;q<dim;q++)
// 								vec2[p-1] += ( nodes[ map_in->npmap[h][n] * dim + q ] - stdvec[dmb[0]*dim+q] )
// 										* ( stdvec[dmb[p]*dim+q] - stdvec[dmb[0]*dim+q] );
// 							sum = 0.0;
// 							for(q=0;q<dim;q++)
// 								sum += powf( stdvec[dmb[p]*dim+q] - stdvec[dmb[0]*dim+q], 2.0 );
// 							vec2[p-1] /= sqrt( sum );
// 						}

int ltgmap_gen_ltgmap( ltgmap_t *map_in )
{
	long int i,j,k,m,n,nface,idx,jdx,pdx,cur = 0,dim = map_in->basis->dim;
	long int *tmp,cmb[dim],*local = NULL,*localk = NULL;
	tmp = (long int*) malloc( map_in->mesh->nsimp * map_in->basis->pdim * sizeof(long int) );
	if( tmp == NULL )
		return -1;
	else
		map_in->map = tmp;

	/* Iterate through each plane */
	for(i=0;i<map_in->nplane;i++)
	{
		/* Make sure local is big enough */
		if( map_in->psize[i] > cur )
		{
			tmp = (long int*) realloc( local, map_in->psize[i] * sizeof(long int) );
			if( tmp == NULL )
				return -1;
			local = tmp;
			tmp = (long int*) realloc( localk, map_in->psize[i] * sizeof(long int) );
			if( tmp == NULL )
				return -1;
			localk = tmp;
			cur = map_in->psize[i];
		}

		/* Iterate through each simplex containing it */
		for(j=0;j<map_in->mbsize[i];j++)
		{
			/* Search for map_in->plane[i] in each simplex to get local indices */
			if( ltgmap_simplex_local_indices( map_in, map_in->member[i][j], map_in->psize[i], map_in->plane[i], local ) == -1 )
				return -2;

			/* Now calculate combinadic index of local and which local nodes lie in the local plane */
			bubble_sort_long( map_in->psize[i], local );
			idx = combinadic_index( map_in->psize[i], local );
			for(k=1;k<map_in->psize[i];k++)
				idx += binomial( dim + 1, k );

			/* Now find all other local definitions of this plane local to simplices in member and check them too */
			for(k=0;k<j;k++)
			{
				/* Do the same as above for member k */
				if( ltgmap_simplex_local_indices( map_in, map_in->member[i][k], map_in->psize[i], map_in->plane[i], localk ) == -1 )
					return -2;

				/* Calculate the offset as above for k */
				bubble_sort_long( map_in->psize[i], localk );
				jdx = combinadic_index( map_in->psize[i], localk );
				for(m=1;m<map_in->psize[i];m++)
					jdx += binomial( dim + 1, m );

				/* Now use map_in->npeqvs and plane indices idx and jdx */
				/* Calculate the pair combinadic index for the plane pair idx and jdx */
				if( idx < jdx )
					localk[0] = idx, localk[1] = jdx;
				else
					localk[0] = jdx, localk[1] = idx;
				pdx = combinadic_index( 2, localk );
				for(m=1;m<map_in->psize[i];m++)
					pdx += binomial( binomial( dim + 1, m ), 2l );

				/* Now check pair index pdx in map_in->npeqvs for which nodes are equivalent */
				for(m=0;m<map_in->eqsize[pdx];m++)
				{
					/* Look at map_in->npeqvs[pdx]: The list of pairs of
					   equivalent nodes in this pair of planes forms an
					   equivalence class of global nodes which should of
					   course only be given a single global index. */
					map_in->npeqvs[pdx][2*m];
					/* If two local nodes are listed in npeqvs[pdx], then
					   (map_in->member[i][j],npeqvs[pdx][0]) and
					   (map_in->member[i][k],npeqvs[pdx][1]) are equivalent */
				}
			}
		}
	}
}

/* Look at pairs of surfaces of the same dimension */
for(i=1,oset=0,poset=0;i<=dim;oset+=binomial(dim+1,i),poset+=binomial(binomial(dim+1,i),2l),i++)
{
	bin = binomial( dim + 1, i );
	combinadic_init( dim + 1, i, cmb );
	for(j=0;j<bin;j++)
	{
		// The local plane index of cmb
		g = oset + j;

		// Initialize second combinadic
		combinadic_init( dim + 1, i, dmb );

		for(k=0;k<j;k++)
		{
			// The local plane index of dmb
			h = oset + k;

			/* Now make combinations of nodes in each plane */
			for(m=0;m<map_in->npsize[g];m++)
			{
				for(p=0;p<i;p++)
					for(q=0;q<dim;q++)
						smpmat[(dim-p)*dim+q] = stdvec[cmb[p]*dim+q];
				for(p=0,r=0;p<dim+1&&r<dim+1;p++)
				{
					if( simple_find_long( p, i, cmb ) == -1 )
					{
						for(q=0;q<dim;q++)
							smpmat[r*dim+q] = stdvec[p*dim+q];
						++r;
					}
				}
				basis_barycentric( dim, nodes + map_in->npmap[g][m] * dim, smpmat, vec1 );

				for(n=0;n<map_in->npsize[h];n++)
				{
					/* CALCULATE THE BARYCENTRIC COORDINATES */
					for(p=0;p<i;p++)
						for(q=0;q<dim;q++)
							smpmat[(dim-p)*dim+q] = stdvec[dmb[p]*dim+q];
					for(p=0,r=0;p<dim+1&&r<dim+1;p++)
					{
						if( simple_find_long( p, i, dmb ) == -1 )
						{
							for(q=0;q<dim;q++)
								smpmat[r*dim+q] = stdvec[p*dim+q];
							++r;
						}
					}
					basis_barycentric( dim, nodes + map_in->npmap[h][n] * dim, smpmat, vec2 );

					/* Now compare vec2 to vec1 for equality */
					q = 0;
					for(p=0;p<dim+1;p++) // Vector is of dimension i - 1, not i
					{
						if( fabs( vec2[p] - vec1[p] ) > LTGMAP_ZERO_THRESHOLD )
						{
							q = 1;
							break;
						}
					}

					/* If q == 0, then equal, so map (g,m) to (h,n) as equal */
					if( q == 0 )
					{
						pvc[0] = k; pvc[1] = j; // Because k < j
						p = poset + combinadic_index( 2, pvc ); // The correctly offset plane pair combinadic index
						++map_in->eqsize[p];
						map_in->npeqvs[p] = (long int*) realloc( map_in->npeqvs[p], 2 * map_in->eqsize[p] * sizeof(long int) ); // npeqvs[p] is the map for plane pair p, but node indices are still global
						map_in->npeqvs[p][2*map_in->eqsize[p]-2] = map_in->npmap[g][m]; // From plane j
						map_in->npeqvs[p][2*map_in->eqsize[p]-1] = map_in->npmap[h][n]; // From plane k
					}
				}
			}
			combinadic_next( dim + 1, i, dmb );
		}
		combinadic_next( dim + 1, i, cmb );
	}
}

int ltgmap_gen_ltgmap( ltgmap_t *map_in )
{
	long int i,j,k,m,n,p,q,r,s,coset,cidx1,cidx2,fidx1,fidx2,eqx,ngn,cur = 0,dim = map_in->basis->dim,pdim = map_in->basis->pdim;
	long int *tmp,cmb[dim],local[dim],locals[dim],localf[dim];
	int mode,match;

	tmp = (long int*) malloc( map_in->mesh->nsimp * map_in->basis->pdim * sizeof(long int) );
	if( tmp == NULL )
		return -1;
	map_in->map = tmp;
	for(i=0;i<map_in->mesh->nsimp*pdim;i++)
		map_in->map[i] = -2; // Initial value is -2 since -1 means something else (as yet unassigned, but needed)

	/* The number of global nodes */
	ngn = -1;

	/* Iterate through each plane */
	for(i=0;i<map_in->nplane;i++)
	{
		/* Form the combinadic offset for size map_in->psize[i] */
		coset = 0;
		for(j=1;j<map_in->psize[i];j++)
			coset += binomial( dim + 1, j );

		/* Iterate through each simplex containing plane i */
		for(j=0;j<map_in->mbsize[i];j++)
		{
			/* Get the definition of plane i in simplex map_in->member[i][j] */
			if( ltgmap_simplex_local_indices( map_in, map_in->member[i][j], map_in->psize[i], map_in->plane[i], local ) == -1 )
				return -3;
			for(k=0;k<map_in->psize[i];k++)
				locals[k] = local[k];
			bubble_sort_long( map_in->psize[i], locals );
			for(k=0;k<map_in->psize[i];k++)
			{
				m = simple_find_long( locals[k], map_in->psize[i], local );
				if( m == -1 )
					return -4;
				localf[m] = k;
			}

			/* Now calculate factoradic index from localf */
			fidx1 = factoradic_index_long( map_in->psize[i], localf );
			cidx1 = combinadic_index( map_in->psize[i], locals );

			/* In this round definitely assign a global index to node p */
			for(p=0;p<map_in->npsize[coset+cidx1];p++)
			{
				if( map_in->map[map_in->member[i][j]*pdim+map_in->npmap[coset+cidx1][p]] != -2 )
					continue; /* Continue because already found all equivalencies for this one */

				/* If it isn't defined, then define the fucker and take a step forward */
				map_in->map[map_in->member[i][j]*pdim+map_in->npmap[coset+cidx1][p]] = ngn++;

				/* Now start the second round */
				for(k=0;k<map_in->mbsize[i];k++) /* FIGURE OUT IF YOU NEED TO GO PAST j HERE!!! PROBABLY NEED TO!!! */
				{
					if( ltgmap_simplex_local_indices( map_in, map_in->member[i][k], map_in->psize[i], map_in->plane[i], local ) == -1 )
						return -5;
					for(m=0;m<map_in->psize[i];m++)
						locals[m] = local[m];
					bubble_sort_long( map_in->psize[i], locals );
					for(m=0;m<map_in->psize[i];m++)
					{
						n = simple_find_long( locals[m], map_in->psize[i], local );
						if( n == -1 )
							return -6;
						localf[n] = m;
					}
	
					/* Do it again... */
					fidx2 = factoradic_index_long( map_in->psize[i], localf );
					cidx2 = combinadic_index( map_in->psize[i], locals );
	
					/* Now calculate the entry in map_in->npeqvs */
					eqx = ltgmap_face_pair_index( dim, map_in->psize[i], cidx1, cidx2, fidx1, fidx2 );
	
					/* Now look through map_in->npeqvs[eqx] which is of length map_in->eqsize[eqx] for equivalencies */
					for(q=0;q<map_in->npsize[coset+cidx2];q++)
					{
						for(r=0;r<map_in->eqsize[eqx];r++)
						{
							if( map_in->npeqvs[eqx][2*r+0] == map_in->npmap[coset+cidx1][p]
							 && map_in->npeqvs[eqx][2*r+1] == map_in->npmap[coset+cidx2][q] )
							{
								/* Then [coset+cidx2][q] is equivalent to [coset+cidx1][p] */
								map_in->map[map_in->member[i][k]*pdim+map_in->npmap[coset+cidx2][q]]
									= map_in->map[map_in->member[i][j]*pdim+map_in->npmap[coset+cidx1][p]];
								break;
							}
							else
							if( map_in->npeqvs[eqx][2*r+0] == map_in->npmap[coset+cidx2][q]
							 && map_in->npeqvs[eqx][2*r+1] == map_in->npmap[coset+cidx1][p] )
							{
								map_in->map[map_in->member[i][j]*pdim+map_in->npmap[coset+cidx1][p]]
									= map_in->map[map_in->member[i][k]*pdim+map_in->npmap[coset+cidx2][q]];
								break;
							}
						}
					}
				}
			}
		}
	}
	printf( "ngn = %d\n", ngn );
	return 0;
}

/**
 * Special note here: User specifies the DIMENSION here, i.e. mesh_in->dim is NOT
 * the dimension used. However, if dimension less than mesh_in->dim is specified,
 * the first dim_in coordinates are those used. Otherwise write another fucking
 * function. Shit.
 */
double __mesh_simplex_determinant( mesh_t *mesh_in, long int dim_in, long int *vids_in )
{
	dmatrix_t mat;
	long int i,j;
	double prd = 1.0;
	int dim = dim_in;
	mat.start = (double*) malloc( ( dim + 1 ) * ( dim + 1 ) * sizeof(double) );
	mat.rsize = dim + 1;
	mat.nrows = dim + 1;
	mat.ncols = dim + 1;

	for(i=0;i<dim+1;i++)
	{
		for(j=0;j<dim;j++)
			mat.start[j*mat.rsize+i] = mesh_in->vertx[vids_in[i]*dim+j];
	}
	for(j=0;j<dim+1;j++)
		mat.start[dim*mat.rsize+j] = 1.0;

	dmatrix_lu_dcmp( &mat );
	for(i=0;i<dim+1;i++)
		prd *= mat.start[i*mat.rsize+i];

	free( mat.start );

	return prd;
}

/* FIXME: THIS ONE SHOULD BE REPLACED TOO */
double mesh_vert_above_face_vert( mesh_t *mesh_in, double *vert_in, long int *vids_in )
{
	/* Equivalent to plane_test as written previously */
	long int i,j,vidx;
	double sum;
	dmatrix_t mat;
	dvector_t vec;
	mat.start = (double*) malloc( ( mesh_in->dim + 1 ) * ( mesh_in->dim + 1 ) * sizeof(double) );
	mat.rsize = mesh_in->dim + 1;
	mat.nrows = mesh_in->dim + 1;
	mat.ncols = mesh_in->dim + 1;
	vec.start = (double*) malloc( ( mesh_in->dim + 1 ) * sizeof(double) );
	vec.rsize = 1;
	vec.length = mesh_in->dim + 1;

	for(i=0;i<mesh_in->dim+1;i++)
	{
		vidx = vids_in[i];
		mat.start[i*mat.rsize] = 1.0;
		for(j=0;j<mesh_in->dim;j++)
			mat.start[i*mat.rsize+j+1] = mesh_in->vertx[ mesh_in->dim * vidx + j ];
	}
	for(i=0;i<mesh_in->dim;i++)
		vec.start[i*vec.rsize] = 1.0;
	vec.start[mesh_in->dim*vec.rsize] = 0.0;
	dmatrix_lu_dcmp( &mat );
	dmatrix_fb_sub_vrhs( &mat, &vec );

	/* Now calculate the interpolant's value at vert_in */
	sum = vec.start[0];
	for(i=0;i<mesh_in->dim;i++)
		sum += vec.start[(i+1)*vec.rsize] * vert_in[i];

	free( mat.start );
	free( vec.start );

	return sum;
}

/**
 * This function tests whether an input node is above a plane of specified
 * index by generating (on-the-fly) a linear interpolation function which is
 * zero at the simplex's vertex which is not part of the dim nodes comprising
 * the face, and unity at all vertices in the plane. This means that if the
 * node in vert_in is above the face, mesh_vert_above_face will return
 * a value greater than 1.0. Otherwise, the value will be less than 1.0.
 * \param mesh_in The mesh object on which to operate
 * \param vert_in The vertex to test
 * \param fid_in This is the index of the unfinished face in ufaces
 */
/* FIXME: THIS FUNCTION SHOULD BE REPLACED BY pt_surf_test in dlinalg.h */
double mesh_vert_above_face( mesh_t *mesh_in, double *vert_in, long int fid_in )
{
	/* Equivalent to plane_test as written previously */
	long int i;
	double sum;
	long int vidx[mesh_in->dim+1];
	for(i=0;i<mesh_in->dim+1;i++)
		vidx[i] = mesh_in->ufaces[fid_in*(mesh_in->dim+1)+i];
	sum = mesh_vert_above_face_vert( mesh_in, vert_in, vidx );
	return sum;
}

/**
 * This function returns 1 if point_in lies to the RIGHT of face_in, which
 * is defined by dim_in points.
 */
int dmatrix_pt_surf_test( long int dim_in, double *face_in, double *point_in )
{
	long int i,j;
	double det1,det2;
	dmatrix_t mat1,mat2;

	mat1.start = (double*) malloc( dim_in * dim_in * sizeof(double) );
	mat1.rsize = dim_in;
	mat1.nrows = dim_in;
	mat1.ncols = dim_in;
	mat2.start = (double*) malloc( ( dim_in + 1 ) * ( dim_in + 1 ) * sizeof(double) );
	mat2.rsize = dim_in + 1;
	mat2.nrows = dim_in + 1;
	mat2.ncols = dim_in + 1;

	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			mat1.start[j*mat1.rsize+i] = face_in[i*dim_in+j];
	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			mat2.start[j*mat2.rsize+i] = face_in[i*dim_in+j];
	for(j=0;j<dim_in;j++)
		mat2.start[j*mat2.rsize+dim_in] = point_in[j];
	for(i=0;i<dim_in+1;i++)
		mat2.start[dim_in*mat2.rsize+i] = 1.0;

	dmatrix_lu_dcmp( &mat1 );
	det1 = 1.0;
	for(j=0;j<dim_in;j++)
		det1 *= mat1.start[j*mat1.rsize+j];
	det2 = 1.0;
	dmatrix_lu_dcmp( &mat2 );
	for(j=0;j<dim_in+1;j++)
		det2 *= mat2.start[j*mat2.rsize+j];
	free( mat1.start );
	free( mat2.start );
	if( det1 < 0.0 )
	{
		if( det2 > 0.0 )
			return 1;
		else
			return 0;
	}
	else
	{
		if( det2 < 0.0 )
			return 1;
		else
			return 0;
	}
	return -1; // Should never get here; this is a bad bad error
}

double dmatrix_simplex_determ( long int dim_in, double *face_in, double *point_in )
{
	long int i,j;
	double det2;
	dmatrix_t mat2;

	mat2.start = (double*) malloc( ( dim_in + 1 ) * ( dim_in + 1 ) * sizeof(double) );
	mat2.rsize = dim_in + 1;
	mat2.nrows = dim_in + 1;
	mat2.ncols = dim_in + 1;

	for(i=0;i<dim_in;i++)
		for(j=0;j<dim_in;j++)
			mat2.start[j*mat2.rsize+i] = face_in[i*dim_in+j];
	for(j=0;j<dim_in;j++)
		mat2.start[j*mat2.rsize+dim_in] = point_in[j];
	for(i=0;i<dim_in+1;i++)
		mat2.start[dim_in*mat2.rsize+i] = 1.0;

	det2 = 1.0;
	dmatrix_lu_dcmp( &mat2 );
	for(j=0;j<dim_in+1;j++)
		det2 *= mat2.start[j*mat2.rsize+j];

	free( mat2.start );
	return det2;
}

int mesh_gen_del_wrap( mesh_t *mesh_in )
{
	long int h,i,j,k,m,n,p;
	long int vids[mesh_in->dim+1];
	long int wids[mesh_in->dim+1];
	long int xids[mesh_in->dim+1];
	long int vec[mesh_in->dim+1];
	int ck1,ck2;
	double surf,circ;

	for(h=0;mesh_in->nunff>0;h++)
	{
		/* Remove a face from mesh_in->ufaces, f */
		for(i=0;i<mesh_in->dim+1;i++)
			vids[i] = mesh_in->ufaces[i];
		mesh_del_unf_face( mesh_in, 0 );

		/* Now find vertex which finishes f */
		for(i=0;i<mesh_in->nvert;i++)
		{
			/* Make sure node i is not already in vids[0:dim] */
			ck1 = 0;
			for(j=0;j<mesh_in->dim+1;j++)
				if( vids[j] == i )
			{
				ck1 = 1;
				break;
			}
			if( ck1 == 1 )
				continue; // i is already in the set, so continue
			surf = mesh_vert_above_face_vert( mesh_in, mesh_in->vertx + i*mesh_in->dim, vids );
			if( surf > 1.0 )
			{
				for(k=0;k<mesh_in->dim;k++)
					wids[k] = vids[k];
				wids[mesh_in->dim] = i;
				m = 0;
				for(j=0;j<mesh_in->nvert;j++)
				{
					ck1 = 0;
					for(k=0;k<mesh_in->dim+1;k++)
						if( wids[k] == j )
					{
						ck1 = 1;
						break;
					}
					if( ck1 == 1 )
						continue;
					circ = mesh_vert_in_circsph_vert( mesh_in, mesh_in->vertx + j*mesh_in->dim, wids );
					if( circ < 0.0 )
					{
						m = 1;
						break;
					}
				}
				if( m == 0 ) // => no node in its circumsphere
				{
					/* Then the circumsphere is empty, so add to simplices */
					p = mesh_simpx_index( mesh_in, wids );
					if( p == -1 )
						mesh_app_simpx( mesh_in, wids );
					else
						continue;

					/* Now check the faces other than f */
					combinadic_init( mesh_in->dim + 1, mesh_in->dim, vec );
					for(j=0,k=mesh_in->dim-1;j<binomial(mesh_in->dim+1,mesh_in->dim)-1;j++,k--)
					{
						combinadic_next( mesh_in->dim + 1, mesh_in->dim, vec );
						for(m=0;m<mesh_in->dim;m++)
							xids[m] = wids[vec[m]];
						xids[mesh_in->dim] = wids[k];

						/* Check xids against ufaces and add accordingly */
						p = mesh_ufaces_index( mesh_in, xids );
						if( p == -1 )
							mesh_app_unf_face( mesh_in, xids );
						else
							mesh_del_unf_face( mesh_in, p );
					}
					break;
				}
			}
		}
	}

	return 0;
}

int list_init( list_t *list_in )
{
	long int *tmp;

	list_in->start = NULL;
	list_in->pointer = NULL;
	list_in->size = 0;
	list_in->alloc = LIST_SIZE_INC;
	tmp = (long int*) realloc( list_in->pointer, LIST_SIZE_INC * sizeof(long int) );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	return 0;
}

int list_free( list_t *list_in )
{
	if( list_in->start != NULL )
		free( list_in->pointer );
	else
		return -1;
	return 0;
}

int list_reset( list_t *list_in )
{
	if( list_in->size < 0 )
		return -1;
	memmove( list_in->pointer, list_in->start, list_in->size * sizeof(long int) );
	list_in->start = list_in->pointer;
	return 0;
}

int list_shift( list_t *list_in, int shift_in )
{
	if( list_in->size <= 0
	 || ( ( list_in->start + shift_in ) < list_in->pointer )
	 || ( ( list_in->start + shift_in ) > ( list_in->pointer + list_in->size ) ) )
		return -1; /* Meaning shift takes pointer out of range */
	if( shift_in == 0 )
		return 0;
	list_in->start = (long int*) memmove( list_in->start + shift_in, list_in->start, list_in->size * sizeof(long int) );
	return 0;
}

/**
 * NOTE: list_t::start should always be the pointer to the first non-null
 * entry in the allocated memory at list_t::pointer. When more memory is
 * needed and list_grow is called, then the data starting at list_t::start
 * is copied back to the beginning in list_t::pointer. Also, list_grow is
 * only called if there is not enough memory available.
 * @param list_in The list to use
 * @param minc_in The number of memory units by which to increment
 */
int list_grow( list_t *list_in, int minc_in )
{
	long int *tmp;
	list_reset( list_in );
	tmp = (long int*) realloc( list_in->pointer, ( list_in->alloc + minc_in ) * sizeof(long int) );
	if( tmp == NULL )
		return -1;
	list_in->pointer = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc += minc_in;

	return list_in->alloc;
}

int list_shrink( list_t *list_in, int mdec_in )
{
	long int *tmp;
	list_reset( list_in ); /* Force list_in->start to coincide with list_in->pointer */
	tmp = (long int*) realloc( list_in->pointer, ( list_in->alloc - mdec_in ) * sizeof(long int) );
	if( tmp == NULL )
		return -1;
	list_in->start = tmp;
	list_in->start = list_in->pointer;
	list_in->alloc -= mdec_in;

	return list_in->alloc;
}

int list_prepend( list_t *list_in, long int val_in )
{
	if( list_in->size + 1 > list_in->alloc )
	{
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
		list_shift( list_in, 1 );
	}
	else
		if( list_in->start <= list_in->pointer ) /* Meaning there is no memory allocated to the left, so make some */
			list_shift( list_in, 1 );

	--list_in->start;
	list_in->start[0] = val_in;
	++list_in->size;

	return list_in->size;
}

int list_append( list_t *list_in, long int val_in )
{
	/* If there isn't enough memory allocated */
	if( list_in->size + 1 > list_in->alloc )
		if( list_grow( list_in, LIST_SIZE_INC ) == -1 )
			return -1;
	else /* If there is enough, then either there is some at the beginning or some at the end */
		if( list_in->start > list_in->pointer ) /* If there is some at the beginning */
			list_shift( list_in, 1 ); /* Shift to the right by sizeof(long) bytes */

	list_in->start[list_in->size] = val_in;
	++list_in->size;

	return list_in->size;
}

long int list_popf( list_t *list_in )
{
	long int out;
	if( list_in->size <= 0 )
		return LIST_NULL_PTR;
	else
		out = list_in->start[0];
	if( list_in->size > 1 )
		++list_in->start;
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		list_shrink( list_in, LIST_SIZE_INC );
	return out;
}

long int list_popb( list_t *list_in )
{
	long int out;

	if( list_in->size <= 0 )
		return LIST_NULL_PTR;
	else
		out = list_in->start[list_in->size-1];
	--list_in->size;
	if( list_in->alloc - list_in->size >= LIST_SIZE_INC )
		list_shrink( list_in, LIST_SIZE_INC );
	return out;
}
