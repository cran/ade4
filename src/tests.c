#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

double betweenvar (double **x, double *w, double *indica);
void reduction (double **tab, double **resu, double *pl, double *pc, int *rang1);
double inerbetweennopondco (double *pl, int moda, double *indica, double **tab);
double inerbetween (double *pl, double *pc, int moda, double *indica, double **tab);
void vecpermut (double *A, int *num, double *B);
void matpermut (double **A, int *num, double **B);
double alea (void);
void aleapermutvec (double *a);
void aleapermutmat (double **a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);

/*********************************************/
void testmantel(int *npermut1,
				int *lig1,
				double *init11,
				double *init21,
				
				double *inersim)
{
/* Declarations de variables C locales */

	int			i, j, k, lig, i0, j0, npermut, *numero, isel;
	double		**m1, **m2;
	double		trace, trace0, moy1, moy2, car1, car2, a0;

/* Allocation memoire pour les variables C locales */

	npermut = *npermut1;
	lig = *lig1;

	taballoc(&m1, lig, lig);
	taballoc(&m2, lig, lig);
	vecintalloc (&numero, lig);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=lig; i++) {
		for (j=1; j<=lig; j++) {
			m1[i][j] = init11[k];
			k = k + 1;
		}
	}

	k = 0;
	for (i=1; i<=lig; i++) {
		for (j=1; j<=lig; j++) {
			m2[i][j] = init21[k];
			k = k + 1;
		}
	}

	trace=0;
	moy1 = 0; moy2=0; car1 = 0; car2 = 0;
	for (i=1; i<=lig; i++) {
		for (j=1; j<=lig; j++) {
			trace = trace + m1[i][j]*m2[i][j];
			if (j>i) {
				moy1 = moy1 + m1[i][j];
				moy2 = moy2 + m2[i][j];
				car1 = car1 + m1[i][j]*m1[i][j];
				car2 = car2 + m2[i][j]*m2[i][j];
			}
		}
	}
	trace = trace/2;
	a0 = trace - moy1*moy2*2/lig/(lig-1);
	a0 = a0/ sqrt ( (double) (car1 - moy1*moy1*2/lig/(lig-1)) );
	a0 = a0/ sqrt ( (double) (car2 - moy2*moy2*2/lig/(lig-1)) );
	trace = a0;
	
	inersim[0] = a0;

	for (isel=1; isel<=npermut; isel++) {
		getpermutation (numero, isel);
		trace0=0;
		for (i=1; i<=lig; i++) {
			i0 = numero[i];
			for (j=1; j<=lig; j++) {
				j0 = numero[j];
				trace0 = trace0 + m1[i][j]*m2[i0][j0];
			}
		}
		trace0 = trace0/2;
		a0 = trace0 - moy1*moy2*2/lig/(lig-1);
		a0 = a0/ sqrt ( (double) (car1 - moy1*moy1*2/lig/(lig-1)) );
		a0 = a0/ sqrt ( (double) (car2 - moy2*moy2*2/lig/(lig-1)) );
		inersim[isel] = a0;
	}

	freetab(m1);
	freetab(m2);
	freeintvec(numero);	
}

/*********************************************/
void testprocuste(	int *npermut1,
				int *lig1,
				int *c11,
				int *c21,
				double *init11,
				double *init21,
				
				double *inersim)
{
/* Declarations de variables C locales */

	int			i, j, k, res, lig, c1, c2, npermut, rang, *numero;
	double		**tabperm, **init1, **init2, tinit, tsim;
	double		**cov, **w, *valpro, *tvecsim;

/* Allocation memoire pour les variables C locales */

	npermut = *npermut1;
	lig = *lig1;
	c1 = *c11;
	c2 = *c21;

	if (c1<=c2) {
		taballoc(&tabperm, lig, c1);
		taballoc(&init1, lig, c1);
		taballoc(&init2, lig, c2);
	} else {
		taballoc(&tabperm, lig, c2);
		taballoc(&init1, lig, c2);
		taballoc(&init2, lig, c1);

		res=c1;
		c1=c2;
		c2=res;
	}	
	taballoc(&cov, c1, c2);
	taballoc(&w, c1, c1);
	vecalloc(&valpro,c1);
	vecintalloc (&numero, lig);
	vecalloc(&tvecsim, npermut);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=lig; i++) {
		for (j=1; j<=c1; j++) {
			init1[i][j] = init11[k];
			k = k + 1;
		}
	}

	k = 0;
	for (i=1; i<=lig; i++) {
		for (j=1; j<=c2; j++) {
			init2[i][j] = init21[k];
			k = k + 1;
		}
	}

/* Calculs */

	tinit = 0;
	prodmatAtBC (init1, init2, cov);
	prodmatAAtB (cov,w);
	DiagobgComp(c1, w, valpro, &rang);
	for (i=1;i<=rang;i++) {
		tinit=tinit+sqrt(valpro[i]);
	}
	
	for (k=1; k<=npermut; k++) {
	
		getpermutation (numero,k);
		matpermut (init1, numero, tabperm);
		
		prodmatAtBC (tabperm, init2, cov);
		prodmatAAtB (cov,w);
		DiagobgComp(c1, w, valpro, &rang);
		tsim=0;
		for (i=1;i<=rang;i++) {
			tsim=tsim+sqrt(valpro[i]);
		}
		tvecsim[k] = tsim;
	}

	inersim[0] = tinit;

	for (k=1; k<=npermut; k++) {
		inersim[k] = tvecsim[k];
	}	
	
	freetab(tabperm);
	freetab(cov);
	freetab(init1);
	freetab(init2);
	freetab(w);
	freevec(tvecsim);
	freevec(valpro);
	freeintvec(numero);
}

/*********************************************/
void testdiscrimin(	int *npermut,
				double *rank,
				double *pl1,
				int *npl,
				int *moda1,
				double *indica1,
				int *nindica,
				double *tab1, int *il1, int *ic1,
				double *inersim)
{
/* Declarations de variables C locales */

	int			l1, c1;
	
	double	**tab, **tabp, *pl, *plp, *indica, rang;
	int		moda, i, j, k, *numero;

/* Allocation memoire pour les variables C locales */

	l1 = *il1;
	c1 = *ic1;
	moda = *moda1;
	rang = *rank;

	vecalloc (&pl, *npl);
	vecalloc (&plp, *npl);
	vecalloc (&indica, *nindica);
	taballoc (&tab, l1, c1);
	taballoc (&tabp, l1, c1);
	vecintalloc(&numero, l1);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			tab[i][j] = tab1[k];
			k = k + 1;
		}
	}
	for (i=1; i<=*npl; i++) {
		pl[i] = pl1[i-1];
	}
	for (i=1; i<=*nindica; i++) {
		indica[i] = indica1[i-1];
	}

/* Calculs
	inertie initiale est stockee dans le premier element du vecteur
	des simulations */

	inersim[0] = betweenvar(tab, pl, indica)/rang;

	for (k=1; k<=*npermut; k++) {
		getpermutation (numero, k);
		matpermut (tab, numero, tabp);
		vecpermut (pl, numero, plp);
		inersim[k] = betweenvar (tabp, plp, indica)/rang;
	}

	freevec (pl);
	freevec (plp);
	freevec (indica);
	freetab (tab);
	freetab (tabp);
	freeintvec (numero);
}

/*********************************************/
double betweenvar (double **tab, double *pl, double *indica)
{
	double	*m, s, bvar, *indicaw;
	int		i, j, l1, c1, ncla, icla;
	
	l1 = tab[0][0];
	c1 = tab[1][0];

	ncla = indica[1];
	for (i=1;i<=l1;i++) {
		if (indica[i] > ncla) ncla = indica[i];
	}
	
	vecalloc(&m, ncla);
	vecalloc(&indicaw, ncla);
	
	bvar = 0;
	for (j=1;j<=c1;j++) {

		for (i=1;i<=ncla;i++) {
			m[i] = 0;
			indicaw[i] = 0;
		}

		for (i=1;i<=l1;i++) {
			icla = indica[i];
			indicaw[icla] = indicaw[icla] + pl[i];
			m[icla] = m[icla] + tab[i][j] * pl[i];
		}
		
		s = 0;
		for (i=1;i<=ncla;i++) {
			s = s + m[i] * m[i] / indicaw[i];
		}
		
		bvar = bvar + s;
	}
	
	freevec(m);
	freevec(indicaw);

	return (bvar);
}

/************************************/
void testinter(	int *npermut,
				double *pl1,
				int *npl,
				double *pc1,
				int *npc,
				int *moda1,
				double *indica1,
				int *nindica,
				double *tab1, int *l1, int *c1,
				double *inersim)
{
/* Declarations de variables C locales */

	double	**tab, **tabp, *pl, *plp, *pc, *indica;
	int		moda, i, j, k;
	int		*numero;
	
/* Allocation memoire pour les variables C locales */

	moda = *moda1;
	vecalloc (&pl, *npl);
	vecalloc (&plp, *npl);
	vecalloc (&pc, *npc);
	vecalloc (&indica, *nindica);
	taballoc (&tab, *l1, *c1);
	taballoc (&tabp, *l1, *c1);
	vecintalloc(&numero, *l1);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=*l1; i++) {
		for (j=1; j<=*c1; j++) {
			tab[i][j] = tab1[k];
			k = k + 1;
		}
	}
	for (i=1; i<=*npl; i++) {
		pl[i] = pl1[i-1];
	}
	for (i=1; i<=*npc; i++) {
		pc[i] = pc1[i-1];
	}
	for (i=1; i<=*nindica; i++) {
		indica[i] = indica1[i-1];
	}

/* Calculs
	inertie initiale est stockee dans le premier element du vecteur
	des simulations */
	
	inersim[0] = inerbetween (pl, pc, moda, indica, tab);

	for (k=1; k<=*npermut; k++) {
		getpermutation (numero,k);
		matpermut (tab, numero, tabp);
		vecpermut (pl, numero, plp);
		inersim[k] = inerbetween (plp, pc, moda, indica, tabp);
	}

	freetab(tab);
	freetab(tabp);
	freevec(pl);
	freevec(plp);
	freevec(pc);
	freevec(indica);
	freeintvec(numero);
}

/************************************/
double inerbetween (double *pl, double *pc, int moda, double *indica, double **tab)
{
	int i, j, k, l1, rang;
	double poi,  inerb, a0, a1, s1;
	double **moy;
	double *pcla;
	
	l1 = tab[0][0];
	rang = tab[1][0];
	taballoc (&moy, moda, rang);
	vecalloc (&pcla, moda);
	
	for (i=1;i<=l1;i++) { 
		k = (int) indica[i];
		poi = pl[i];
		pcla[k]=pcla[k]+poi;
	}

	
	for (i=1;i<=l1;i++) {
		k = (int) indica[i];
		poi = pl[i];
		for (j=1;j<=rang;j++) {
			moy[k][j] = moy[k][j] + tab[i][j]*poi;
		}
	}
	
	for (k=1;k<=moda;k++) { 
		a0 = pcla[k];
		for (j=1;j<=rang;j++) {
			moy[k][j] = moy[k][j]/a0;
		}
	}

	inerb = 0;
	for (i=1;i<=moda;i++) {
		a1 = pcla[i];
		for (j=1;j<=rang;j++) {
			s1 = moy[i][j];
			inerb = inerb + s1 * s1 *a1 * pc[j];
		}
	}
	freetab (moy);
	freevec (pcla);
	return inerb;
	
}

/*****************/
void testertrace (	int *npermut,
				double *pc1r, int *npc1,
				double *pc2r, int *npc2,
				double *tab1r, int *l1r, int *c1r,
				double *tab2r, int *l1r1, int *c2r,
				double *inersimul)
{

/* Declarations des variables C locales */

	double	**X1, **X2, *pc1, *pc2, **cov;
	int		i, j, k, l1, c1, c2;
	double	poi, inertot, s1, inersim;
	int		*numero;
	
/* On recopie les objets R dans les variables C locales */

	l1 = *l1r;
	c1 = *c1r;
	c2 = *c2r;
		
/* Allocation memoire pour les variables C locales */

	vecalloc (&pc1, *npc1);
	vecalloc (&pc2, *npc2);
	vecintalloc(&numero, l1);
	taballoc (&X1, l1, c1);
	taballoc (&X2, l1, c2);
	taballoc(&cov, c2, c1);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			X1[i][j] = tab1r[k];
			k = k + 1;
		}
	}
	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c2; j++) {
			X2[i][j] = tab2r[k];
			k = k + 1;
		}
	}
	for (i=1; i<=*npc1; i++) {
		pc1[i] = pc1r[i-1];
	}
	for (i=1; i<=*npc2; i++) {
		pc2[i] = pc2r[i-1];
	}

/* Calculs */

	for (j=1;j<=c1;j++) {
		poi = sqrt(pc1[j]);
		for (i=1; i<=l1;i++) {
			X1[i][j]=X1[i][j]*poi;
		}
	}
	for (j=1;j<=c2;j++) {
		poi = sqrt(pc2[j]);
		for (i=1; i<=l1;i++) {
			X2[i][j]=X2[i][j]*poi;
		}
	}
	
	prodmatAtBC (X2, X1, cov);
	
	inertot = 0;
	for (i=1;i<=c2;i++) {
		for (j=1;j<=c1;j++) {
			s1 = cov[i][j];
			inertot = inertot + s1 * s1;
		}
	}
	inertot = inertot / l1 / l1;
	inersimul[0] = inertot;
	
	for (k=1; k<=*npermut; k++) {
		getpermutation (numero,k);
		prodmatAtBrandomC (X2, X1, cov, numero);

		inersim = 0;
		for (i=1;i<=c2;i++) {
			for (j=1;j<=c1;j++) {
				s1 = cov[i][j];
				inersim = inersim + s1 * s1;
			}
		}
		inersimul[k] = inersim / l1 / l1;
	}
	
	freevec (pc1);
	freevec (pc2);
	freeintvec (numero);
	freetab (X1);
	freetab (X2);
	freetab (cov);
}

/*****************/
void testertracenu (	int *npermut,
				double *pc1r, int *npc1,
				double *pc2r, int *npc2,
				double *plr, int *npl,
				double *tab1r, int *l1r, int *c1r,
				double *tab2r, int *l1r1, int *c2r,
				double *tabinit1r,
				double *tabinit2r,
				int *typ1r,
				int *typ2r,
				double *inersimul)
{
/* Declarations des variables C locales */

	double	**X1, **X2, **init1, **init2, *pc1, *pc2, *pl, **cov;
	int		i, j, k, l1, c1, c2;
	double	poi, inertot, s1, inersim, a1;
	int		*numero1, *numero2;
	char	typ1[3], typ2[3];
	
/* On recopie les objets R dans les variables C locales */

	l1 = *l1r;
	c1 = *c1r;
	c2 = *c2r;
	strncpy(typ1, (char const *) *typ1r, 2);
	strncpy(typ2, (char const *) *typ2r, 2);
	typ1[2] = 0;
	typ2[2] = 0;

/* Allocation memoire pour les variables C locales */

	vecalloc (&pc1, *npc1);
	vecalloc (&pc2, *npc2);
	vecalloc (&pl, l1);
	vecintalloc (&numero1, l1);
	vecintalloc (&numero2, l1);
	taballoc (&X1, l1, c1);
	taballoc (&X2, l1, c2);
	taballoc (&init1, l1, c1);
	taballoc (&init2, l1, c2);
	taballoc (&cov, c2, c1);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			init1[i][j] = tab1r[k];
			k = k + 1;
		}
	}
	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c2; j++) {
			init2[i][j] = tab2r[k];
			k = k + 1;
		}
	}
	for (i=1; i<=*npc1; i++) {
		pc1[i] = pc1r[i-1];
	}
	for (i=1; i<=*npc2; i++) {
		pc2[i] = pc2r[i-1];
	}
	for (i=1; i<=*npl; i++) {
		pl[i] = plr[i-1];
	}
	
/* Calculs */

	inertot = 0;
	for (i=1; i<=l1;i++) {
		poi = pl[i];
		for (j=1;j<=c1;j++) {
			init1[i][j]=init1[i][j]*poi;
		}
	}

	prodmatAtBC (init2, init1, cov);
	
	for (i=1;i<=c2;i++) {
		a1 = pc2[i];
		for (j=1;j<=c1;j++) {
			s1 = cov[i][j];
			inertot = inertot + s1 * s1 * a1 * pc1[j];
		}
	}

	inersimul[0] = inertot;

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			init1[i][j] = tabinit1r[k];
			k = k + 1;
		}
	}
	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c2; j++) {
			init2[i][j] = tabinit2r[k];
			k = k + 1;
		}
	}
	
	for (k=1; k<=*npermut; k++) {
	
		getpermutation (numero1,k);
		getpermutation (numero2,2*k);
		
		matpermut (init1, numero1, X1);
		matpermut (init2, numero2, X2);
		
		matcentrage (X1, pl, typ1);
		matcentrage (X2, pl, typ2);

		for (i=1; i<=l1;i++) {
			poi = pl[i];
			for (j=1;j<=c1;j++) {
				X1[i][j]=X1[i][j]*poi;
			}
		}

		prodmatAtBC (X2, X1, cov);

		inersim = 0;
		for (i=1;i<=c2;i++) {
			a1 = pc2[i];
			for (j=1;j<=c1;j++) {
				s1 = cov[i][j];
				inersim = inersim + s1 * s1 * a1 * pc1[j];
			}
		}
		inersimul[k] = inersim;
	}
	freevec (pc1);
	freevec (pc2);
	freevec (pl);
	freeintvec (numero1);
	freeintvec (numero2);
	freetab (X1);
	freetab (X2);
	freetab (init1);
	freetab (init2);
	freetab (cov);
}


/*****************/
void testertracenubis (	int *npermut,
				double *pc1r, int *npc1,
				double *pc2r, int *npc2,
				double *plr, int *npl,
				double *tab1r, int *l1r, int *c1r,
				double *tab2r, int *l1r1, int *c2r,
				double *tabinit1r,
				double *tabinit2r,
				int *typ1r,
				int *typ2r,
				int *ntabr,
				double *inersimul)

{
/* Declarations des variables C locales */

	double	**X1, **X2, **init1, **init2, *pc1, *pc2, *pl, **cov;
	int		i, j, k, l1, c1, c2;
	double	poi, inertot, s1, inersim, a1;
	int		*numero1, *numero2, ntab;
	char	typ1[3], typ2[3];

/* On recopie les objets R dans les variables C locales */

	l1 = *l1r;
	c1 = *c1r;
	c2 = *c2r;
	ntab = *ntabr;
	strncpy(typ1, (char const *) *typ1r, 2);
	strncpy(typ2, (char const *) *typ2r, 2);
	typ1[2] = 0;
	typ2[2] = 0;
			
/* Allocation memoire pour les variables C locales */

	vecalloc (&pc1, *npc1);
	vecalloc (&pc2, *npc2);
	vecalloc (&pl, l1);
	vecintalloc (&numero1, l1);
	vecintalloc (&numero2, l1);
	taballoc (&X1, l1, c1);
	taballoc (&X2, l1, c2);
	taballoc (&init1, l1, c1);
	taballoc (&init2, l1, c2);
	taballoc (&cov, c2, c1);

/* On recopie les objets R dans les variables C locales */

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			init1[i][j] = tab1r[k];
			k = k + 1;
		}
	}
	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c2; j++) {
			init2[i][j] = tab2r[k];
			k = k + 1;
		}
	}
	for (i=1; i<=*npc1; i++) {
		pc1[i] = pc1r[i-1];
	}
	for (i=1; i<=*npc2; i++) {
		pc2[i] = pc2r[i-1];
	}
	for (i=1; i<=*npl; i++) {
		pl[i] = plr[i-1];
	}
	
	inertot = 0;
	for (i=1; i<=l1;i++) {
		poi = pl[i];
		for (j=1;j<=c1;j++) {
			init1[i][j]=init1[i][j]*poi;
		}
	}

	prodmatAtBC (init2, init1, cov);
	
	for (i=1;i<=c2;i++) {
		a1 = pc2[i];
		for (j=1;j<=c1;j++) {
			s1 = cov[i][j];
			inertot = inertot + s1 * s1 * a1 * pc1[j];
		}
	}
	
	inersimul[0] = inertot;

	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c1; j++) {
			X1[i][j] = tab1r[k];
			k = k + 1;
		}
	}
	for (i=1; i<=l1;i++) {
		poi = pl[i];
		for (j=1;j<=c1;j++) {
			X1[i][j]=X1[i][j]*poi;
		}
	}
	
	k = 0;
	for (i=1; i<=l1; i++) {
		for (j=1; j<=c2; j++) {
			X2[i][j] = tab2r[k];
			k = k + 1;
		}
	}
	for (i=1; i<=l1;i++) {
		poi = pl[i];
		for (j=1;j<=c2;j++) {
			X2[i][j]=X2[i][j]*poi;
		}
	}
	
	if (ntab == 1) {
		k = 0;
		for (i=1; i<=l1; i++) {
			for (j=1; j<=c2; j++) {
				init2[i][j] = tabinit2r[k];
				k = k + 1;
			}
		}
	} else {
		k = 0;
		for (i=1; i<=l1; i++) {
			for (j=1; j<=c1; j++) {
				init1[i][j] = tabinit1r[k];
				k = k + 1;
			}
		}
	}

	for (k=1; k<=*npermut; k++) {
	
		if (ntab == 1) {
			getpermutation (numero2,k);
			matpermut (init2, numero2, X2);
			matcentrage (X2, pl, typ2);
		} else {
			getpermutation (numero1,k);
			matpermut (init1, numero1, X1);
			matcentrage (X1, pl, typ1);
		}
		
		prodmatAtBC (X2, X1, cov);

		inersim = 0;
		for (i=1;i<=c2;i++) {
			a1 = pc2[i];
			for (j=1;j<=c1;j++) {
				s1 = cov[i][j];
				inersim = inersim + s1 * s1 * a1 * pc1[j];
			}
		}
		inersimul[k] = inersim;
	}
	freevec (pc1);
	freevec (pc2);
	freevec (pl);
	freeintvec (numero1);
	freeintvec (numero2);
	freetab (X1);
	freetab (X2);
	freetab (init1);
	freetab (init2);
	freetab (cov);
}

/**************************/
double alea (void)
{
	double w;
	w = ((double) rand())/ (double)RAND_MAX;
	return (w);
}

/*************************/
void aleapermutvec (double *a)
{
	// permute au hasard les ÚlÚments du vecteur a
	// Manly p. 42 Le vecteur est modifiÚ
	// from Knuth 1981 p. 139
	int lig, i,j, k;
	double z;
	
	lig = a[0];
	for (i=1; i<=lig-1; i++) {
		j=lig-i+1;
		k = (int) (j*alea()+1);
		//k = (int) (j*genrand()+1);
		if (k>j) k=j;
		z = a[j];
		a[j]=a[k];
		a[k] = z;
	}
}

/*************************/
void aleapermutmat (double **a)
{
	// permute au hasard les lignes du tableau a
	// Manly p. 42 le tableau est modifiÚ
	int lig, i,j, col, n, k;
	double z;

	lig = a[0][0];
	col = a[1][0];
	for (i=1; i<=lig-1; i++) {
		j=lig-i+1;
		k = (int) (j*alea ()+1);
		//k = (int) (j*genrand()+1);
		if (k>j) k=j;
		for (n=1; n<=col; n++) {
			z = a[j][n];
			a[j][n]=a[k][n];
			a[k][n] = z;
		}
	}
}

/*******************/	
void matpermut (double **A, int *num, double **B)
{
/*---------------------------------------
* A est une matrice n-p
* B est une matrice n-p
* num est une permutation alÚatoire des n premiers entiers
* B contient en sortie les lignes de A permutÚes
* ---------------------------------------*/

	int lig, col,lig1, col1, lig2, i, j, k;
	
	lig = A[0][0];
	col = A[1][0];
	lig1 = B[0][0];
	col1 = B[1][0];
	lig2 = num[0];
	
	
	if ( (lig!=lig1) || (col!=col1) || (lig!=lig2) ) {
		return;
	}
	
	for (i=1; i<=lig; i++) {
		k=num[i];
		for (j=1; j<=col; j++) {
			B[i][j] = A[k][j];
		}
	}
}

/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
* A est un vecteur n elements
* B est une vecteur n elements
* num est une permutation alŽatoire des n premiers entiers
* B contient en sortie les elements de A permutŽes
* ---------------------------------------*/

	int lig, lig1, lig2, i, k;
	
	lig = A[0];
	lig1 = B[0];
	lig2 = num[0];
	
	
	if ( (lig!=lig1) || (lig!=lig2) ) {
		//err_message ("Illegal parameters (vecpermut)");
		//closelisting();
	}
	
	for (i=1; i<=lig; i++) {
		k=num[i];
		B[i] = A[k];
	}
}

/*******************/	
void matcentrage (double **A, double *poili, char *typ)
{
	
	if (strcmp (typ,"nc") == 0) {
		return;
	} else if (strcmp (typ,"cm") == 0) {
		matmodifcm (A, poili);
		return;
	} else if (strcmp (typ,"cn") == 0) {
		matmodifcn (A, poili);
		return;
	} else if (strcmp (typ,"cp") == 0) {
		matmodifcp (A, poili);
		return;
	} else if (strcmp (typ,"cs") == 0) {
		matmodifcs (A, poili);
		return;
	} else if (strcmp (typ,"fc") == 0) {
		matmodiffc (A, poili);
		return;
	} else if (strcmp (typ,"fl") == 0) {
		matmodifcm (A, poili);
		return;
	}
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, m colonnes
* disjonctif complet
* poili est un vecteur n composantes
* la procedure retourne tab centre par colonne 
* pour la ponderation poili (somme=1)
* centrage type correspondances multiples
--------------------------------------------------*/
{
	double		poid;
	int 			i, j, l1, m1;
	double		*poimoda;
	double		x, z;

	l1 = tab[0][0];
	m1 = tab[1][0];
	vecalloc(&poimoda, m1);


	for (i=1;i<=l1;i++) {
		poid = poili[i];
		for (j=1;j<=m1;j++) {
			poimoda[j] = poimoda[j] + tab[i][j] * poid;
		}
	}
	
	for (j=1;j<=m1;j++) {
		x = poimoda[j];
		if (x==0) {
			for (i=1;i<=l1;i++) tab[i][j] = 0;
		} else {
		
			for (i=1;i<=l1;i++) {
				z = tab[i][j]/x - 1.0;
				tab[i][j] = z;
			}
		}
	}
	freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab norme par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
	double		poid, x, z, y, v2;
	int 			i, j, l1, c1;
	double		*moy, *var;

	l1 = tab[0][0];
	c1 = tab[1][0];

	vecalloc(&moy, c1);
	vecalloc(&var, c1);


/*--------------------------------------------------
* calcul du tableau centre/norme
--------------------------------------------------*/

	for (i=1;i<=l1;i++) {
		poid = poili[i];
		for (j=1;j<=c1;j++) {
			moy[j] = moy[j] + tab[i][j] * poid;
		}
	}
	
	for (i=1;i<=l1;i++) {
		poid=poili[i];
		for (j=1;j<=c1;j++) {
			x = tab[i][j] - moy[j];
			var[j] = var[j] + poid * x * x;
		}
	}
	
	for (j=1;j<=c1;j++) {
		v2 = var[j];
		if (v2<=0) v2 = 1;
		v2 = sqrt(v2);
		var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
		x = moy[i];
		y = var[i];
		for (j=1;j<=l1;j++) {
			z = tab[j][i] - x;
			z = z / y;
			tab[j][i] = z;
		}
	}
	
	freevec(moy);
	freevec(var);
	
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab standardise par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
	double		poid, x, z, y, v2;
	int 			i, j, l1, c1;
	double		*moy, *var;

	l1 = tab[0][0];
	c1 = tab[1][0];

	vecalloc(&var, c1);


/*--------------------------------------------------
* calcul du tableau standardise
--------------------------------------------------*/

	for (i=1;i<=l1;i++) {
		poid=poili[i];
		for (j=1;j<=c1;j++) {
			x = tab[i][j];
			var[j] = var[j] + poid * x * x;
		}
	}
	
	for (j=1;j<=c1;j++) {
		v2 = var[j];
		if (v2<=0) v2 = 1;
		v2 = sqrt(v2);
		var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
		x = moy[i];
		y = var[i];
		for (j=1;j<=l1;j++) {
			z = tab[j][i];
			z = z / y;
			tab[j][i] = z;
		}
	}
	freevec(var);
}

/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, p colonnes
* poili est un vecteur n composantes
* la procedure retourne tab centre par colonne 
* pour la ponderation poili (somme=1)
--------------------------------------------------*/
{
	double		poid;
	int 			i, j, l1, c1;
	double		*moy, x, z;

	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&moy, c1);


/*--------------------------------------------------
* calcul du tableau centre
--------------------------------------------------*/

	for (i=1;i<=l1;i++) {
		poid = poili[i];
		for (j=1;j<=c1;j++) {
			moy[j] = moy[j] + tab[i][j] * poid;
		}
	}
	
	
	for (i=1;i<=c1;i++) {
		x = moy[i];
		for (j=1;j<=l1;j++) {
			z = tab[j][i] - x;
			tab[j][i] = z;
		}
	}
	freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
* tab est un tableau n lignes, m colonnes
* de nombres positifs ou nuls
* poili est un vecteur n composantes
* la procedure retourne tab centre doublement 
* pour la ponderation poili (somme=1)
* centrage type correspondances simples
--------------------------------------------------*/
{
	double		poid;
	int 			i, j, l1, m1;
	double		*poimoda;
	double		x, z;

	l1 = tab[0][0];
	m1 = tab[1][0];
	vecalloc(&poimoda, m1);


	for (i=1;i<=l1;i++) {
		x = 0;
		for (j=1;j<=m1;j++) {
			x = x + tab[i][j];
		}
		if (x!=0) {
			for (j=1;j<=m1;j++) {
				tab[i][j] = tab[i][j]/x;
			}
		}	
	}

	for (i=1;i<=l1;i++) {
		poid = poili[i];
		for (j=1;j<=m1;j++) {
			poimoda[j] = poimoda[j] + tab[i][j] * poid;
		}
	}
	
	for (j=1;j<=m1;j++) {
		x = poimoda[j];
		if (x==0) {
			//err_message("column has a nul weight (matmodiffc)");
		}
		
		for (i=1;i<=l1;i++) {
			z = tab[i][j]/x - 1.0;
			tab[i][j] = z;
		}
	}
	freevec (poimoda);
}

/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
* affectation d'une permutation alÚatoire des n premiers entiers 
* dans dans un vecteur d'entiers de dimension n
* vecintalloc prÚalable exigÚ
* *numero est un vecteur d'entier
* repet est un entier qui peut prendre une valeur arbitraire
* utilise dans le germe du generateur de nb pseudo-aleatoires
* si on l'incremente dans des appels repetes (e.g. simulation) garantit
* que deux appels donnent deux resultats distincts (seed=clock+repet)
------------------------*/
{
	int i, n, seed;
	int *alea;
	
	n=numero[0];
	vecintalloc (&alea,n);
	
	/*-------------
	* numerotation dans numero
	-----------*/
	for (i=1;i<=n;i++) {
		numero[i]=i;
	}
	
	/*-------------
	* affectation de nombres aleatoires dans alea
	----------------*/
	seed = clock();
	seed = seed + repet;
	srand(seed);
	for (i=1;i<=n;i++) {
		alea[i]=rand();
	}
	
	trirapideint (alea , numero, 1, n);
	freeintvec (alea);
}

/*****************************************/
void trirapideint (int *x , int *num, int gauche, int droite)
{
	int j, dernier, milieu, t;
	
	if ( (droite-gauche)<=0) return;
	
	milieu = (gauche+droite)/2;
	trirapideintswap (x, gauche, milieu);
	trirapideintswap (num, gauche, milieu);
	
	t=x[gauche];
	dernier=gauche;
	for (j = gauche+1; j<=droite; j++) {
		if (x[j] < t) {
			dernier = dernier + 1;
			trirapideintswap (x, dernier, j);	
			trirapideintswap (num, dernier, j);
		}
	}
	trirapideintswap (x, gauche, dernier);
	trirapideintswap (num, gauche, dernier);
	
	trirapideint (x, num, gauche, dernier-1);
	trirapideint (x, num, dernier+1, droite);
		
}

/**************************************/
void trirapideintswap (int *v, int i, int j)
{
	int provi;
	
	provi=v[i];
	v[i]=v[j];
	v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
* Racine carree des elements d'un vecteur
--------------------------------------------------*/
{
	int i, c1;
	double v2;
	
	c1 = v1[0];
	
	for (i=1;i<=c1;i++) {
		v2 = v1[i];
		// if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)");
		v2 = sqrt(v2);
		v1[i] = v2;
	}
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
* Diagonalisation
* T. FOUCART Analyse factorielle de tableaux multiples,
* Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
* de LEBART et coll.
--------------------------------------------------*/
{
	double			*s;
	double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
	double			dble;
	int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
	
	vecalloc(&s, n0);
	a = 0.000000001;
	ni = 100;
	if (n0 == 1) {
		d[1] = w[1][1];
		w[1][1] = 1.0;
		*rang = 1;
		freevec (s);
		return;
	}
	
	for (i2=2;i2<=n0;i2++) {
		
		b=0.0;
		c=0.0;
		i=n0-i2+2;
		k=i-1;
		if (k < 2) goto Et1;
		for (l=1;l<=k;l++) {
			c = c + fabs((double) w[i][l]);
		}
		if (c != 0.0) goto Et2;
		
Et1:	s[i] = w[i][k];
		goto Etc;
		
Et2:	for (l=1;l<=k;l++) {
			x = w[i][l] / c;
			w[i][l] = x;
			b = b + x * x;
		}
		xp = w[i][k];
		ix = 1;
		if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
		dble = b;
		dble = -sqrt(dble);
		q = dble * ix;

		s[i] = c * q;
		b = b - xp * q;
		w[i][k] = xp - q;
		xp = 0;
		for (m=1;m<=k;m++) {
			w[m][i] = w[i][m] / b / c;
			q = 0;
			for (l=1;l<=m;l++) {
				q = q + w[m][l] * w[i][l];
			}
			m1 = m + 1;
			if (k < m1) goto Et3;
			for (l=m1;l<=k;l++) {
				q = q + w[l][m] * w[i][l];
			}
			
Et3:		s[m] = q / b;
			xp = xp + s[m] * w[i][m];
		}
		bp = xp * 0.5 / b;
		for (m=1;m<=k;m++) {
			xp = w[i][m];
			q = s[m] - bp * xp;
			s[m] = q;
			for (l=1;l<=m;l++) {
				w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
			}
		}
		for (l=1;l<=k;l++) {
			w[i][l] = c * w[i][l];
		}
		
Etc:	d[i] = b;
	} /* for (i2=2;i2<n0;i2++) */
	
	s[1] = 0.0;
	d[1] = 0.0;
	
	for (i=1;i<=n0;i++) {
		
		k = i - 1;
		if (d[i] == 0.0) goto Et4;
		for (m=1;m<=k;m++) {
			q = 0.0;
			for (l=1;l<=k;l++) {
				q = q + w[i][l] * w[l][m];
			}
			for (l=1;l<=k;l++) {
				w[l][m] = w[l][m] - q * w[l][i];
			}
		}
		
Et4:	d[i] = w[i][i];
		w[i][i] = 1.0;
		if (k < 1) goto Et5;
		for (m=1;m<=k;m++) {
			w[i][m] = 0.0;
			w[m][i] = 0.0;
		}

Et5:;
	}
	
	for (i=2;i<=n0;i++) {
		s[i-1] = s[i];
	}
	s[n0] = 0.0;
	
	for (k=1;k<=n0;k++) {

		m = 0;

Et6: 	for (j=k;j<=n0;j++) {
			if (j == n0) goto Et7;
			ab = fabs((double) s[j]);
			ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
			if (ab < ep) goto Et7;
		}
	
Et7: 	isnou = 1;
		h = d[k];
		if (j == k) goto Eta;
		if (m < ni) goto Etd;
		
		//err_message("Error: can't compute matrix eigenvalues");
		
Etd:	m = m + 1;
		q = (d[k+1]-h) * 0.5 / s[k];
		
/*		t = sqrt(q * q + 1.0); */
		dble = q * q + 1.0;
		dble = sqrt(dble);
		t = dble;
		
		if (q < 0.0) isnou = -1;
		q = d[j] - h + s[k] / (q + t * isnou);
		u = 1.0;
		v = 1.0;
		h = 0.0;
		jk = j-k;
		for (ijk=1;ijk<=jk;ijk++) {
			i = j - ijk;
			xp = u * s[i];
			b = v * s[i];
			if (fabs((double) xp) < fabs((double) q)) goto Et8;
			u = xp / q;
			
/*			t = sqrt(u * u + 1); */
			dble = u * u + 1.0;
			dble = sqrt(dble);
			t = dble;
			
			s[i+1] = q * t;
			v = 1 / t;
			u = u * v;
			goto Et9;

Et8:		v = q / xp;

/*			t = sqrt(1 + v * v); */
			dble = 1.0 + v * v;
			dble = sqrt(dble);
			t = dble;
			
			s[i+1] = t * xp;
			u = 1 / t;
			v = v * u;

Et9:
			q = d[i+1] - h;
			t = (d[i] - q) * u + 2.0 * v * b;
			h = u * t;
			d[i+1] = q + h;
			q = v * t - b;
			for (l=1;l<=n0;l++) {
				xp = w[l][i+1];
				w[l][i+1] = u * w[l][i] + v * xp;
				w[l][i] = v * w[l][i] - u * xp;
			}
		}
		d[k] = d[k] - h;
		s[k] = q;
		s[j] = 0.0;
		
		goto Et6;

Eta:;
	} /* for (k=1;k<=n0;k++) */
	
	for (ij=2;ij<=n0;ij++) {
		
		i = ij - 1;
		l = i;
		h = d[i];
		for (m=ij;m<=n0;m++) {
			if (d[m] >= h) {
				l = m;
				h = d[m];
			}
		}
		if (l == i) {
			goto Etb;
		} else {
			d[l] = d[i];
			d[i] = h;
		}
		for (m=1;m<=n0;m++) {
			h = w[m][i];
			w[m][i] = w[m][l];
			w[m][l] = h;
		}

Etb:;
	} /* for (ij=2;ij<=n0;ij++) */

//final:;
	*rang = 0;
	for (i=1;i<=n0;i++) {
		/*
		if (d[i] / d[1] < 0.00001) d[i] = 0.0;
		if (d[i] != 0.0) *rang = *rang + 1;
		*/
		if (d[i] > 0.0) *rang = *rang + 1;
	}
	freevec(s);
} /* DiagoCompbg */

/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Produit matriciel AB
--------------------------------------------------*/
{
	int j, k, i, lig, col, col2;
	double s;
	
	lig = a[0][0];
	col = a[1][0];
	
	col2 = b[1][0];

	for (i=1;i<=lig;i++) {
		for (k=1;k<=col2;k++) {
			s = 0;
			for (j=1;j<=col;j++) {
				s = s + a[i][j] * b[j][k];
			}
		c[i][k] = s;
		}		
	}
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Produit matriciel AtA
--------------------------------------------------*/
{
	int j, k, i, lig, col;
	double s;
	
	lig = a[0][0];
	col = a[1][0];

	for (j=1;j<=col;j++) {
		for (k=j;k<=col;k++) {
			s = 0;
			for (i=1;i<=lig;i++) {
				s = s + a[i][k] * a[i][j];
			}
		b[j][k] = s;
		b[k][j] = s;
		}		
	}
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
* Produit matriciel AtB
--------------------------------------------------*/
{
	int j, k, i, lig, col, col2;
	double s;
	
	lig = a[0][0];
	col = a[1][0];
	
	col2 = b[1][0];

	for (j=1;j<=col;j++) {
		for (k=1;k<=col2;k++) {
			s = 0;
			for (i=1;i<=lig;i++) {
				s = s + a[i][j] * b[i][k];
			}
		c[j][k] = s;
		}		
	}
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
* Produit matriciel B = AAt
--------------------------------------------------*/
{
	int j, k, i, lig, col;
	double s;
	
	lig = a[0][0];
	col = a[1][0];

	for (j=1;j<=lig;j++) {
		for (k=j;k<=lig;k++) {
			s = 0;
			for (i=1;i<=col;i++) {
				s = s + a[j][i] * a[k][i];
			}
		b[j][k] = s;
		b[k][j] = s;
		}		
	}
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
* Produit matriciel AtB
* les lignes de B sont permutÚes par la permutation permut
--------------------------------------------------*/
{
	int j, k, i, i0, lig, col, col2;
	double s;
	
	lig = a[0][0];
	col = a[1][0];
	
	col2 = b[1][0];

	for (j=1;j<=col;j++) {
		for (k=1;k<=col2;k++) {
			s = 0;
			for (i=1;i<=lig;i++) {
				i0 = permut[i];
				s = s + a[i][j] * b[i0][k];
			}
		c[j][k] = s;
		}		
	}
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau (l1, c1)
--------------------------------------------------*/
{
	int i, j;
	
	if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
		for (i=0;i<=l1;i++) {
			if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
				return;
				for (j=0;j<i;j++) {
					free(*(*tab+j));
				}
			}
		}
	}

	**(*tab) = l1;
	**(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
* Allocation de memoire pour un vecteur de longueur n
--------------------------------------------------*/
{
	if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
		**vec = n;
		return;
	} else {
		return;
	}
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
* Allocation de memoire pour un vecteur d'entiers de longueur n
--------------------------------------------------*/
{
	if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
		**vec = n;
		return;
	} else {
		return;
	}
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
* Allocation de memoire dynamique pour un tableau (l1, c1)
--------------------------------------------------*/
{
	int 	i, n;
	
	n = *(*(tab));
	for (i=0;i<=n;i++) {
			free((char *) *(tab+i) );
	}
	free((char *) tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
* liberation de memoire pour un vecteur
--------------------------------------------------*/
{
	free((char *) vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* liberation de memoire pour un vecteur
--------------------------------------------------*/
{
	
	free((char *) vec);
	
}


