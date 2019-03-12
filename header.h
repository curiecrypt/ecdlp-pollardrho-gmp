#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

typedef struct {
	mpz_t x;
	mpz_t y;
	mpz_t z;
} POINT_t[1], *POINT;
typedef struct {
	POINT_t P;
	POINT_t Q;
	mpz_t A;
	mpz_t p;
	mpz_t k;
} SMUL_t[1], *SMUL;
typedef struct {
	POINT_t Xd;
	POINT_t Xdd;
	POINT_t P;
	POINT_t Q;
	mpz_t cd;
	mpz_t dd;
	mpz_t cdd;
	mpz_t ddd;
	mpz_t p;
	mpz_t n;
	mpz_t k;
} LOOP_t[1], *LOOP;
void sum(SMUL w) {
	mpz_t x1, x2, y1, y2, z1, z2, temp1, temp2, temp3, temp4, tempx1, tempy1;

	//initialize the temporary variables for operations
	mpz_init(temp1);
	mpz_init(temp2);
	mpz_init(temp3);
	mpz_init(temp4);
	mpz_init(tempx1);
	mpz_init(tempy1);

	mpz_init(x1);
	mpz_init(x2);
	mpz_init(y1);
	mpz_init(y2);
	mpz_init(z1);
	mpz_init(z2);

	//setting the variables to multiplicands
	mpz_set(x1, w->P->x);
	mpz_set(x2, w->Q->x);
	mpz_set(y1, w->P->y);
	mpz_set(y2, w->Q->y);
	mpz_set(z1, w->P->z);
	mpz_set(z2, w->Q->z);

	if (mpz_cmp_ui(x1, 0) == 0 && mpz_cmp_ui(y1, 1) == 0) {
		mpz_set(w->P->x, x2);
		mpz_set(w->P->y, y2);
		mpz_set(w->P->z, z2);
	} else if (mpz_cmp_ui(x2, 0) == 0 && mpz_cmp_ui(y2, 1) == 0) {
		mpz_set(w->P->x, x1);
		mpz_set(w->P->y, y1);
		mpz_set(w->P->z, z1);
	} else if (mpz_cmp(x1, x2) == 0) {
		if (mpz_cmp(y1, y2) != 0) {
			mpz_set_ui(w->P->x, 0);
			mpz_set_ui(w->P->y, 1);
			mpz_set_ui(w->P->z, 0);
		} else if (mpz_cmp_ui(y2, 0) == 0) {
			mpz_set_ui(w->P->x, 0);
			mpz_set_ui(w->P->y, 1);
			mpz_set_ui(w->P->z, 0);
		} else
			mpz_mul(temp1, x1, x1); //x^2
		mpz_mod(temp1, temp1, w->p);
		mpz_mul_ui(temp2, temp1, 3); //3*x^2
		mpz_mod(temp2, temp2, w->p);
		mpz_add(temp1, temp2, w->A); //3*x^2+A
		mpz_mod(temp1, temp1, w->p);
		mpz_mul_ui(temp2, y1, 2); //2*y1
		mpz_mod(temp2, temp2, w->p);
		mpz_invert(temp3, temp2, w->p);
		mpz_mod(temp3, temp3, w->p);
		mpz_mul(temp2, temp1, temp3); //(3*x^2+A)/(2*y1)
		mpz_mod(temp2, temp2, w->p);
		mpz_mul(temp1, temp2, temp2); //((3*x^2+A)/(2*y1))^2
		mpz_mod(temp1, temp1, w->p);
		mpz_mul_ui(temp2, x1, 2); //2*x1
		mpz_mod(temp2, temp2, w->p);
		mpz_sub(temp3, temp1, temp2); //((3*x^2+A)/(2*y1))^2-(2*x1)

		mpz_mod(tempx1, temp3, w->p);

		mpz_mul(temp1, x1, x1); //x^2
		mpz_mod(temp1, temp1, w->p);
		mpz_mul_ui(temp2, temp1, 3); //3*x^2
		mpz_mod(temp2, temp2, w->p);
		mpz_add(temp1, temp2, w->A); //3*x^2+A
		mpz_mod(temp1, temp1, w->p);
		mpz_mul_ui(temp2, y1, 2); //2*y1
		mpz_mod(temp2, temp2, w->p);
		mpz_invert(temp3, temp2, w->p);
		mpz_mod(temp3, temp3, w->p);
		mpz_mul(temp2, temp1, temp3); //(3*x^2+A)/(2*y1)
		mpz_mod(temp2, temp2, w->p);

		mpz_sub(temp1, x1, tempx1); //x1-x3
		mpz_mod(temp1, temp1, w->p);

		mpz_mul(temp3, temp2, temp1); //(3*x^2+A)/(2*y1)*(x1-x3)
		mpz_mod(temp3, temp3, w->p);
		mpz_sub(temp2, temp3, y1); //(3*x^2+A)/(2*y1)*(x1-x3)-y1
		mpz_mod(tempy1, temp2, w->p);

		mpz_set(w->P->x, tempx1);
		mpz_set(w->P->y, tempy1);
		mpz_set_ui(w->P->z, 1);
	} else {
		mpz_sub(temp1, y1, y2); //y1-y2
		mpz_mod(temp1, temp1, w->p);
		mpz_sub(temp2, x1, x2); //x1-x2
		mpz_mod(temp2, temp2, w->p);
		mpz_invert(temp3, temp2, w->p);
		mpz_mod(temp3, temp3, w->p);
		mpz_mul(temp2, temp1, temp3); //(y1- y2)/(x1-x2)
		mpz_mod(temp2, temp2, w->p);
		mpz_mul(temp1, temp2, temp2); //((y1- y2)/(x1-x2))^2
		mpz_mod(temp1, temp1, w->p);
		mpz_sub(temp2, temp1, x1); //((y1- y2)/(x1-x2))^2 - x1
		mpz_mod(temp2, temp2, w->p);
		mpz_sub(temp1, temp2, x2); //((y1- y2)/(x1-x2))^2 - x1- x2
		mpz_mod(tempx1, temp1, w->p);

		mpz_sub(temp1, y1, y2); //y1-y2
		mpz_mod(temp1, temp1, w->p);
		mpz_sub(temp2, x1, x2); //x1-x2
		mpz_mod(temp2, temp2, w->p);
		mpz_invert(temp3, temp2, w->p);
		mpz_mod(temp3, temp3, w->p);
		mpz_mul(temp2, temp1, temp3); //(y1- y2)/(x1-x2)
		mpz_mod(temp2, temp2, w->p);
		mpz_sub(temp1, x1, tempx1); //x1-x3
		mpz_mod(temp1, temp1, w->p);
		mpz_mul(temp3, temp2, temp1); //(y1- y2)/(x1-x2) * (x1-x3)
		mpz_mod(temp3, temp3, w->p);
		mpz_sub(temp2, temp3, y1); //(y1- y2)/(x1-x2) * (x1-x3) -y1
		mpz_mod(tempy1, temp2, w->p);

		mpz_set(w->P->x, tempx1);
		mpz_set(w->P->y, tempy1);
		mpz_set_ui(w->P->z, 1);

		mpz_clear(temp1);
		mpz_clear(temp2);
		mpz_clear(temp3);
		mpz_clear(temp4);
		mpz_clear(tempx1);
		mpz_clear(tempy1);

		mpz_clear(x1);
		mpz_clear(x2);
		mpz_clear(y1);
		mpz_clear(y2);
		mpz_clear(z1);
		mpz_clear(z2);

	}
}
void dbl(SMUL wm) {
	mpz_t x1, x2, y1, y2, z1, z2, temp1, temp2, temp3, temp4, tempx1, tempy1,
			tempmod1, tempmod2, tempmod3;

	//initialize the temporary variables for operations
	mpz_init(temp1);
	mpz_init(temp2);
	mpz_init(temp3);
	mpz_init(temp4);
	mpz_init(tempx1);
	mpz_init(tempy1);
	mpz_init(tempmod1);
	mpz_init(tempmod2);
	mpz_init(tempmod3);

	mpz_init(x1);
	mpz_init(x2);
	mpz_init(y1);
	mpz_init(y2);
	mpz_init(z1);
	mpz_init(z2);

	//setting the variables to multiplicands
	mpz_set(x1, wm->P->x);
	mpz_set(x2, wm->P->x);
	mpz_set(y1, wm->P->y);
	mpz_set(y2, wm->P->y);
	mpz_set(z1, wm->P->z);
	mpz_set(z2, wm->P->z);

	if (mpz_cmp_ui(x1, 0) == 0 && mpz_cmp_ui(y1, 1) == 0) {
		mpz_set(wm->P->x, x2);
		mpz_set(wm->P->y, y2);
		mpz_set(wm->P->z, z2);
	} else if (mpz_cmp_ui(x2, 0) == 0 && mpz_cmp_ui(y2, 1) == 0) {
		mpz_set(wm->P->x, x1);
		mpz_set(wm->P->y, y1);
		mpz_set(wm->P->z, z1);
	} else if (mpz_cmp(x1, x2) == 0) {
		if (mpz_cmp(y1, y2) != 0) {
			mpz_set_ui(wm->P->x, 0);
			mpz_set_ui(wm->P->y, 1);
			mpz_set_ui(wm->P->z, 0);
		} else if (mpz_cmp_ui(y2, 0) == 0) {
			mpz_set_ui(wm->P->x, 0);
			mpz_set_ui(wm->P->y, 1);
			mpz_set_ui(wm->P->z, 0);
		} else

			mpz_mul(temp1, x1, x1);
		mpz_mod(tempmod1, temp1, wm->p);
		mpz_mul_ui(temp2, tempmod1, 3);
		mpz_mod(tempmod2, temp2, wm->p);
		mpz_add(temp3, tempmod2, wm->A);
		mpz_mod(tempmod3, temp3, wm->p);
		mpz_mul_ui(temp1, y1, 2);
		mpz_mod(tempmod1, temp1, wm->p);
		mpz_invert(temp2, tempmod1, wm->p);
		mpz_mul(temp3, tempmod3, temp2);
		mpz_mod(tempmod3, temp3, wm->p);

		mpz_mul(temp1, tempmod3, tempmod3);
		mpz_mod(tempmod1, temp1, wm->p);
		mpz_mul_ui(temp2, x1, 2);
		mpz_mod(tempmod2, temp2, wm->p);
		mpz_sub(temp3, tempmod1, tempmod2);
		mpz_mod(tempx1, temp3, wm->p);

		mpz_sub(temp1, x1, tempx1);
		mpz_mod(tempmod1, temp1, wm->p);
		mpz_mul(temp2, tempmod3, tempmod1);
		mpz_mod(tempmod2, temp2, wm->p);
		mpz_sub(temp1, tempmod2, y1);
		mpz_mod(tempy1, temp1, wm->p);

		mpz_set(wm->P->x, tempx1);
		mpz_set(wm->P->y, tempy1);
		mpz_set_ui(wm->P->z, 1);

	} else {
		mpz_sub(temp1, y1, y2); //y1-y2
		mpz_mod(temp1, temp1, wm->p);
		mpz_sub(temp2, x1, x2); //x1-x2
		mpz_mod(temp2, temp2, wm->p);
		mpz_invert(temp3, temp2, wm->p);
		mpz_mod(temp3, temp3, wm->p);
		mpz_mul(temp2, temp1, temp3); //(y1- y2)/(x1-x2)
		mpz_mod(temp2, temp2, wm->p);
		mpz_mul(temp1, temp2, temp2); //((y1- y2)/(x1-x2))^2
		mpz_mod(temp1, temp1, wm->p);
		mpz_sub(temp2, temp1, x1); //((y1- y2)/(x1-x2))^2 - x1
		mpz_mod(temp2, temp2, wm->p);
		mpz_sub(temp1, temp2, x2); //((y1- y2)/(x1-x2))^2 - x1- x2
		mpz_mod(tempx1, temp1, wm->p);

		mpz_sub(temp1, y1, y2); //y1-y2
		mpz_mod(temp1, temp1, wm->p);
		mpz_sub(temp2, x1, x2); //x1-x2
		mpz_mod(temp2, temp2, wm->p);
		mpz_invert(temp3, temp2, wm->p);
		mpz_mod(temp3, temp3, wm->p);
		mpz_mul(temp2, temp1, temp3); //(y1- y2)/(x1-x2)
		mpz_mod(temp2, temp2, wm->p);
		mpz_sub(temp1, x1, tempx1); //x1-x3
		mpz_mod(temp1, temp1, wm->p);
		mpz_mul(temp3, temp2, temp1); //(y1- y2)/(x1-x2) * (x1-x3)
		mpz_mod(temp3, temp3, wm->p);
		mpz_sub(temp2, temp3, y1); //(y1- y2)/(x1-x2) * (x1-x3) -y1
		mpz_mod(tempy1, temp2, wm->p);

		mpz_set(wm->P->x, tempx1);
		mpz_set(wm->P->y, tempy1);
		mpz_set_ui(wm->P->z, 1);

		mpz_clear(temp1);
		mpz_clear(temp2);
		mpz_clear(temp3);
		mpz_clear(temp4);
		mpz_clear(tempx1);
		mpz_clear(tempy1);
		mpz_clear(tempmod1);
		mpz_clear(tempmod2);
		mpz_clear(tempmod3);

		mpz_clear(x1);
		mpz_clear(x2);
		mpz_clear(y1);
		mpz_clear(y2);
		mpz_clear(z1);
		mpz_clear(z2);
	}
}
void mult(SMUL w, mpz_t k) {
	SMUL_t mul;
	mp_bitcnt_t index;
	mp_size_t size;
	mpz_init(mul->P->x);
	mpz_init(mul->P->y);
	mpz_init(mul->P->z);
	mpz_init(mul->Q->x);
	mpz_init(mul->Q->y);
	mpz_init(mul->Q->z);
	mpz_init(mul->p);
	mpz_init(mul->A);

	mpz_set(mul->P->x, w->Q->x);
	mpz_set(mul->P->y, w->Q->y);
	mpz_set(mul->P->z, w->Q->z);
	mpz_set(mul->Q->x, w->P->x);
	mpz_set(mul->Q->y, w->P->y);
	mpz_set(mul->Q->z, w->P->z);
	mpz_set(mul->p, w->p);
	mpz_set(mul->A, w->A);

	size = mpz_sizeinbase(k, 2);
	int i;
	for (i = size - 1; i >= 0; i--) {
		dbl(mul);
		index = mpz_scan1(k, i);
		if (index == i) {
			sum(mul);
		}
	}
	mpz_set(w->Q->x, mul->P->x);
	mpz_set(w->Q->y, mul->P->y);
	mpz_set(w->Q->z, mul->P->z);

	mpz_clear(mul->P->x);
	mpz_clear(mul->P->y);
	mpz_clear(mul->P->z);
	mpz_clear(mul->Q->x);
	mpz_clear(mul->Q->y);
	mpz_clear(mul->Q->z);
	mpz_clear(mul->p);
	mpz_clear(mul->A);
}
unsigned long Hfunc(LOOP v) {
	//return w->P->x->_mp_d[0] % LT_LEN;
	return v->Xd->x->_mp_d[0] & ((unsigned long) 31);
}
unsigned long HfuncD(LOOP v) {
	//return w->P->x->_mp_d[0] % LT_LEN;
	return v->Xdd->x->_mp_d[0] & ((unsigned long) 31);
}
