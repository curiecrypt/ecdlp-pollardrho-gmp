#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "header.h"

#define FP_LEN 1
#define LT_LEN 32 //LOOK UP TABLE LENGTH
#define P_LEN 15 //2^15



int main() {
	SMUL_t w;
	LOOP_t v;

	POINT_t R[LT_LEN];
	mpz_t a[LT_LEN];
	mpz_t b[LT_LEN];
	mpz_t t1, t2;

	POINT_t tempP, tempQ;

	//init workspace w
	mpz_init(w->P->x);
	mpz_init(w->P->y);
	mpz_init(w->P->z);
	mpz_init(w->Q->x);
	mpz_init(w->Q->y);
	mpz_init(w->Q->z);
	mpz_init(w->A);
	mpz_init(w->p);
	mpz_init(w->k);

	//init workspace v
	mpz_init(v->Xd->x);
	mpz_init(v->Xd->y);
	mpz_init(v->Xd->z);
	mpz_init(v->Xdd->x);
	mpz_init(v->Xdd->y);
	mpz_init(v->Xdd->z);
	mpz_init(v->P->x);
	mpz_init(v->P->y);
	mpz_init(v->P->z);
	mpz_init(v->Q->x);
	mpz_init(v->Q->y);
	mpz_init(v->Q->z);
	mpz_init(v->cd);
	mpz_init(v->dd);
	mpz_init(v->cdd);
	mpz_init(v->ddd);
	mpz_init(v->p);
	mpz_init(v->n);
	mpz_init(v->k);

	//init temp. points
	mpz_init(tempP->x);
	mpz_init(tempP->y);
	mpz_init(tempP->z);
	mpz_init(tempQ->x);
	mpz_init(tempQ->y);
	mpz_init(tempQ->z);


	//setting the LOOP workspace 64 bits
	mpz_set_ui(v->P->x, 451396633900755934);
	mpz_set_ui(v->P->y, 1593598289986602382);
	mpz_set_ui(v->P->z, 1);
	mpz_set_ui(v->Q->x, 438274830426242316);
	mpz_set_ui(v->Q->y, 189733902907033708);
	mpz_set_ui(v->Q->z, 1);
	mpz_set_ui(v->p, 2089226132838630061);
	mpz_set_ui(v->n, 2089226130425754689);
	mpz_set_ui(v->k, 1876437903330259560);

	//setting necessary components of SMUL workspace
	mpz_set(w->p, v->p);
	mpz_set(w->k, v->k);
	mpz_set_ui(w->A, 375888055676754941);

	unsigned long int seed;
	gmp_randstate_t r_state;

	seed = 123456789;

	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, seed);

	int i;
	for (i = 0; i < LT_LEN; i++) {
		mpz_init(a[i]);
		mpz_init(b[i]);
		mpz_init(R[i]->x);
		mpz_init(R[i]->y);
		mpz_init(R[i]->z);
		mpz_urandomm(a[i], r_state, v->n);
		mpz_urandomm(b[i], r_state, v->n);

		// R[j] := sum(mult(a[j],P,E,A),mult(b[j],Q,E,A),A,E);********
		mpz_set(w->P->x, v->P->x);
		mpz_set(w->P->y, v->P->y);
		mpz_set(w->P->z, v->P->z);
		mpz_set_ui(w->Q->x, 0);
		mpz_set_ui(w->Q->y, 1);
		mpz_set_ui(w->Q->z, 0);
		mult(w, a[i]);
		mpz_set(tempP->x, w->Q->x);
		mpz_set(tempP->y, w->Q->y);
		mpz_set(tempP->z, w->Q->z);

		mpz_set(w->P->x, v->Q->x);
		mpz_set(w->P->y, v->Q->y);
		mpz_set(w->P->z, v->Q->z);
		mpz_set_ui(w->Q->x, 0);
		mpz_set_ui(w->Q->y, 1);
		mpz_set_ui(w->Q->z, 0);
		mult(w, b[i]);

		mpz_set(w->P->x, tempP->x);
		mpz_set(w->P->y, tempP->y);
		mpz_set(w->P->z, tempP->z);
		sum(w);

		mpz_set(R[i]->x, w->P->x);
		mpz_set(R[i]->y, w->P->y);
		mpz_set(R[i]->z, w->P->z);

		// end**********************************************************
	}
	mpz_urandomm(v->cd, r_state, v->n);
	mpz_urandomm(v->dd, r_state, v->n);

	//Xd := sum(mult(cd,P,E,A), mult(dd,Q,E,A),A,E);******************
	mpz_set(w->P->x, v->P->x);
	mpz_set(w->P->y, v->P->y);
	mpz_set(w->P->z, v->P->z);
	mpz_set_ui(w->Q->x, 0);
	mpz_set_ui(w->Q->y, 1);
	mpz_set_ui(w->Q->z, 0);
	mult(w, v->cd);
	mpz_set(tempP->x, w->Q->x);
	mpz_set(tempP->y, w->Q->y);
	mpz_set(tempP->z, w->Q->z);

	mpz_set(w->P->x, v->Q->x);
	mpz_set(w->P->y, v->Q->y);
	mpz_set(w->P->z, v->Q->z);
	mpz_set_ui(w->Q->x, 0);
	mpz_set_ui(w->Q->y, 1);
	mpz_set_ui(w->Q->z, 0);
	mult(w, v->dd);

	mpz_set(w->P->x, tempP->x);
	mpz_set(w->P->y, tempP->y);
	mpz_set(w->P->z, tempP->z);
	sum(w);

	mpz_set(v->Xd->x, w->P->x);
	mpz_set(v->Xd->y, w->P->y);
	mpz_set(v->Xd->z, w->P->z);

	//end************************************************************

	//Xdd := Xd;cdd := cd; ddd := dd;********************************
	mpz_set(v->Xdd->x, v->Xd->x);
	mpz_set(v->Xdd->y, v->Xd->y);
	mpz_set(v->Xdd->z, v->Xd->z);
	mpz_set(v->cdd, v->cd);
	mpz_set(v->ddd, v->dd);
	//end************************************************************

	//main loop*****************************************************
	unsigned long j;
	do {
		j = Hfunc(v);
		mpz_set(w->P->x, v->Xd->x);
		mpz_set(w->P->y, v->Xd->y);
		mpz_set(w->P->z, v->Xd->z);
		mpz_set(w->Q->x, R[j]->x);
		mpz_set(w->Q->y, R[j]->y);
		mpz_set(w->Q->z, R[j]->z);

		sum(w);
		mpz_set(v->Xd->x, w->P->x);
		mpz_set(v->Xd->y, w->P->y);
		mpz_set(v->Xd->z, w->P->z);

		mpz_add(v->cd, v->cd, a[j]);
		mpz_mod(v->cd, v->cd, v->n);
		mpz_add(v->dd, v->dd, b[j]);
		mpz_mod(v->dd, v->dd, v->n);

		for (i = 1; i <= 2; i++) {
			j = HfuncD(v);

			mpz_set(w->P->x, v->Xdd->x);
			mpz_set(w->P->y, v->Xdd->y);
			mpz_set(w->P->z, v->Xdd->z);
			mpz_set(w->Q->x, R[j]->x);
			mpz_set(w->Q->y, R[j]->y);
			mpz_set(w->Q->z, R[j]->z);
			sum(w);
			mpz_set(v->Xdd->x, w->P->x);
			mpz_set(v->Xdd->y, w->P->y);
			mpz_set(v->Xdd->z, w->P->z);

			mpz_add(v->cdd, v->cdd, a[j]);
			mpz_mod(v->cdd, v->cdd, v->n);
			mpz_add(v->ddd, v->ddd, b[j]);
			mpz_mod(v->ddd, v->ddd, v->n);
		}
	} while (mpz_cmp(v->Xd->x, v->Xdd->x) != 0
			&& mpz_cmp(v->Xd->y, v->Xdd->y) != 0);

	/*mpz_out_str(stdout, 10, v->Xd->x);
	printf("\n");
	mpz_out_str(stdout, 10, v->Xd->y);
	printf("\n");
	mpz_out_str(stdout, 10, v->Xd->z);
	printf("\n");*/

	mpz_init(t1);
	mpz_init(t2);
	if (mpz_cmp(v->dd, v->ddd) == 0) {
		printf("FAILURE !!\n");
	} else {
		mpz_sub(t1, v->cd, v->cdd);
		mpz_sub(t2, v->ddd, v->dd);
		mpz_invert(t2, t2, v->n);
		mpz_mul(t1, t1, t2);
		mpz_mod(t1, t1, v->n);
		printf("Log(P, Q) = ");
		mpz_out_str(stdout, 10, t1);
		printf("\n");
	}

	//Free the occupied memory
	mpz_clear(w->P->x);
	mpz_clear(w->P->y);
	mpz_clear(w->P->z);
	mpz_clear(w->Q->x);
	mpz_clear(w->Q->y);
	mpz_clear(w->Q->z);
	mpz_clear(w->A);
	mpz_clear(w->p);
	mpz_clear(w->k);

	//clear workspace v
	mpz_clear(v->Xd->x);
	mpz_clear(v->Xd->y);
	mpz_clear(v->Xd->z);
	mpz_clear(v->Xdd->x);
	mpz_clear(v->Xdd->y);
	mpz_clear(v->Xdd->z);
	mpz_clear(v->P->x);
	mpz_clear(v->P->y);
	mpz_clear(v->P->z);
	mpz_clear(v->Q->x);
	mpz_clear(v->Q->y);
	mpz_clear(v->Q->z);
	mpz_clear(v->cd);
	mpz_clear(v->dd);
	mpz_clear(v->cdd);
	mpz_clear(v->ddd);
	mpz_clear(v->p);
	mpz_clear(v->n);

	mpz_clear(tempP->x);
	mpz_clear(tempP->y);
	mpz_clear(tempP->z);
	mpz_clear(tempQ->x);
	mpz_clear(tempQ->y);
	mpz_clear(tempQ->z);

	mpz_clear(t1);
	mpz_clear(t2);

	for (j = 0; j < LT_LEN; ++j) {
		mpz_clear(a[j]);
		mpz_clear(b[j]);

		mpz_clear(R[j]->x);
		mpz_clear(R[j]->y);
		mpz_clear(R[j]->z);
	}
}#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "header.h"

#define FP_LEN 1
#define LT_LEN 32 //LOOK UP TABLE LENGTH
#define P_LEN 15 //2^15



int main() {
	SMUL_t w;
	LOOP_t v;

	POINT_t R[LT_LEN];
	mpz_t a[LT_LEN];
	mpz_t b[LT_LEN];
	mpz_t t1, t2;

	POINT_t tempP, tempQ;

	//init workspace w
	mpz_init(w->P->x);
	mpz_init(w->P->y);
	mpz_init(w->P->z);
	mpz_init(w->Q->x);
	mpz_init(w->Q->y);
	mpz_init(w->Q->z);
	mpz_init(w->A);
	mpz_init(w->p);
	mpz_init(w->k);

	//init workspace v
	mpz_init(v->Xd->x);
	mpz_init(v->Xd->y);
	mpz_init(v->Xd->z);
	mpz_init(v->Xdd->x);
	mpz_init(v->Xdd->y);
	mpz_init(v->Xdd->z);
	mpz_init(v->P->x);
	mpz_init(v->P->y);
	mpz_init(v->P->z);
	mpz_init(v->Q->x);
	mpz_init(v->Q->y);
	mpz_init(v->Q->z);
	mpz_init(v->cd);
	mpz_init(v->dd);
	mpz_init(v->cdd);
	mpz_init(v->ddd);
	mpz_init(v->p);
	mpz_init(v->n);
	mpz_init(v->k);

	//init temp. points
	mpz_init(tempP->x);
	mpz_init(tempP->y);
	mpz_init(tempP->z);
	mpz_init(tempQ->x);
	mpz_init(tempQ->y);
	mpz_init(tempQ->z);

	//setting the LOOP workspace 64 bits
	mpz_set_ui(v->P->x, 451396633900755934);
	mpz_set_ui(v->P->y, 1593598289986602382);
	mpz_set_ui(v->P->z, 1);
	mpz_set_ui(v->Q->x, 438274830426242316);
	mpz_set_ui(v->Q->y, 189733902907033708);
	mpz_set_ui(v->Q->z, 1);
	mpz_set_ui(v->p, 2089226132838630061);
	mpz_set_ui(v->n, 2089226130425754689);
	mpz_set_ui(v->k, 1876437903330259560);

	//setting necessary components of SMUL workspace
	mpz_set(w->p, v->p);
	mpz_set(w->k, v->k);
	mpz_set_ui(w->A, 375888055676754941);

	unsigned long int seed;
	gmp_randstate_t r_state;

	seed = 123456789;

	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, seed);

	int i;
	for (i = 0; i < LT_LEN; i++) {
		mpz_init(a[i]);
		mpz_init(b[i]);
		mpz_init(R[i]->x);
		mpz_init(R[i]->y);
		mpz_init(R[i]->z);
		mpz_urandomm(a[i], r_state, v->n);
		mpz_urandomm(b[i], r_state, v->n);

		// R[j] := sum(mult(a[j],P,E,A),mult(b[j],Q,E,A),A,E);********
		mpz_set(w->P->x, v->P->x);
		mpz_set(w->P->y, v->P->y);
		mpz_set(w->P->z, v->P->z);
		mpz_set_ui(w->Q->x, 0);
		mpz_set_ui(w->Q->y, 1);
		mpz_set_ui(w->Q->z, 0);
		mult(w, a[i]);
		mpz_set(tempP->x, w->Q->x);
		mpz_set(tempP->y, w->Q->y);
		mpz_set(tempP->z, w->Q->z);

		mpz_set(w->P->x, v->Q->x);
		mpz_set(w->P->y, v->Q->y);
		mpz_set(w->P->z, v->Q->z);
		mpz_set_ui(w->Q->x, 0);
		mpz_set_ui(w->Q->y, 1);
		mpz_set_ui(w->Q->z, 0);
		mult(w, b[i]);

		mpz_set(w->P->x, tempP->x);
		mpz_set(w->P->y, tempP->y);
		mpz_set(w->P->z, tempP->z);
		sum(w);

		mpz_set(R[i]->x, w->P->x);
		mpz_set(R[i]->y, w->P->y);
		mpz_set(R[i]->z, w->P->z);

		// end**********************************************************
	}
	mpz_urandomm(v->cd, r_state, v->n);
	mpz_urandomm(v->dd, r_state, v->n);

	//Xd := sum(mult(cd,P,E,A), mult(dd,Q,E,A),A,E);******************
	mpz_set(w->P->x, v->P->x);
	mpz_set(w->P->y, v->P->y);
	mpz_set(w->P->z, v->P->z);
	mpz_set_ui(w->Q->x, 0);
	mpz_set_ui(w->Q->y, 1);
	mpz_set_ui(w->Q->z, 0);
	mult(w, v->cd);
	mpz_set(tempP->x, w->Q->x);
	mpz_set(tempP->y, w->Q->y);
	mpz_set(tempP->z, w->Q->z);

	mpz_set(w->P->x, v->Q->x);
	mpz_set(w->P->y, v->Q->y);
	mpz_set(w->P->z, v->Q->z);
	mpz_set_ui(w->Q->x, 0);
	mpz_set_ui(w->Q->y, 1);
	mpz_set_ui(w->Q->z, 0);
	mult(w, v->dd);

	mpz_set(w->P->x, tempP->x);
	mpz_set(w->P->y, tempP->y);
	mpz_set(w->P->z, tempP->z);
	sum(w);

	mpz_set(v->Xd->x, w->P->x);
	mpz_set(v->Xd->y, w->P->y);
	mpz_set(v->Xd->z, w->P->z);

	//end************************************************************

	//Xdd := Xd;cdd := cd; ddd := dd;********************************
	mpz_set(v->Xdd->x, v->Xd->x);
	mpz_set(v->Xdd->y, v->Xd->y);
	mpz_set(v->Xdd->z, v->Xd->z);
	mpz_set(v->cdd, v->cd);
	mpz_set(v->ddd, v->dd);
	//end************************************************************

	//main loop*****************************************************
	unsigned long j;
	do {
		j = Hfunc(v);
		mpz_set(w->P->x, v->Xd->x);
		mpz_set(w->P->y, v->Xd->y);
		mpz_set(w->P->z, v->Xd->z);
		mpz_set(w->Q->x, R[j]->x);
		mpz_set(w->Q->y, R[j]->y);
		mpz_set(w->Q->z, R[j]->z);

		sum(w);
		mpz_set(v->Xd->x, w->P->x);
		mpz_set(v->Xd->y, w->P->y);
		mpz_set(v->Xd->z, w->P->z);

		mpz_add(v->cd, v->cd, a[j]);
		mpz_mod(v->cd, v->cd, v->n);
		mpz_add(v->dd, v->dd, b[j]);
		mpz_mod(v->dd, v->dd, v->n);

		for (i = 1; i <= 2; i++) {
			j = HfuncD(v);

			mpz_set(w->P->x, v->Xdd->x);
			mpz_set(w->P->y, v->Xdd->y);
			mpz_set(w->P->z, v->Xdd->z);
			mpz_set(w->Q->x, R[j]->x);
			mpz_set(w->Q->y, R[j]->y);
			mpz_set(w->Q->z, R[j]->z);
			sum(w);
			mpz_set(v->Xdd->x, w->P->x);
			mpz_set(v->Xdd->y, w->P->y);
			mpz_set(v->Xdd->z, w->P->z);

			mpz_add(v->cdd, v->cdd, a[j]);
			mpz_mod(v->cdd, v->cdd, v->n);
			mpz_add(v->ddd, v->ddd, b[j]);
			mpz_mod(v->ddd, v->ddd, v->n);
		}
	} while (mpz_cmp(v->Xd->x, v->Xdd->x) != 0
			&& mpz_cmp(v->Xd->y, v->Xdd->y) != 0);

	mpz_out_str(stdout, 10, v->Xd->x);
	printf("\n");
	mpz_out_str(stdout, 10, v->Xd->y);
	printf("\n");
	mpz_out_str(stdout, 10, v->Xd->z);
	printf("\n");

	mpz_init(t1);
	mpz_init(t2);
	if (mpz_cmp(v->dd, v->ddd) == 0) {
		printf("FAILURE !!\n");
	} else {
		mpz_sub(t1, v->cd, v->cdd);
		mpz_sub(t2, v->ddd, v->dd);
		mpz_invert(t2, t2, v->n);
		mpz_mul(t1, t1, t2);
		mpz_mod(t1, t1, v->n);
		printf("\n");
		mpz_out_str(stdout, 10, t1);
	}

	//Free the occupied memory
	mpz_clear(w->P->x);
	mpz_clear(w->P->y);
	mpz_clear(w->P->z);
	mpz_clear(w->Q->x);
	mpz_clear(w->Q->y);
	mpz_clear(w->Q->z);
	mpz_clear(w->A);
	mpz_clear(w->p);
	mpz_clear(w->k);

	//clear workspace v
	mpz_clear(v->Xd->x);
	mpz_clear(v->Xd->y);
	mpz_clear(v->Xd->z);
	mpz_clear(v->Xdd->x);
	mpz_clear(v->Xdd->y);
	mpz_clear(v->Xdd->z);
	mpz_clear(v->P->x);
	mpz_clear(v->P->y);
	mpz_clear(v->P->z);
	mpz_clear(v->Q->x);
	mpz_clear(v->Q->y);
	mpz_clear(v->Q->z);
	mpz_clear(v->cd);
	mpz_clear(v->dd);
	mpz_clear(v->cdd);
	mpz_clear(v->ddd);
	mpz_clear(v->p);
	mpz_clear(v->n);

	mpz_clear(tempP->x);
	mpz_clear(tempP->y);
	mpz_clear(tempP->z);
	mpz_clear(tempQ->x);
	mpz_clear(tempQ->y);
	mpz_clear(tempQ->z);

	mpz_clear(t1);
	mpz_clear(t2);

	for (j = 0; j < LT_LEN; ++j) {
		mpz_clear(a[j]);
		mpz_clear(b[j]);

		mpz_clear(R[j]->x);
		mpz_clear(R[j]->y);
		mpz_clear(R[j]->z);
	}
}
