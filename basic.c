#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "flint/fmpz_vec.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpz_poly_mat.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mat.h"
#include "libbgv.h"

double bgv_get_dvn()
{
	return dvn;   
}

void gen_q(fmpz_t q, long len)
{
        mpz_t tmp, hold;
        mpz_init(tmp);
        mpz_init(hold);
        mpz_set_ui(hold, 1);
        mpz_mul_2exp(tmp, hold, len);
        mpz_nextprime(hold, tmp);
        fmpz_set_mpz(q, hold);
        mpz_clear(tmp);
        mpz_clear(hold);
}

void hcrypt_random(mpz_t tmp)
{
	FILE *fp;
	fp = fopen("/dev/urandom", "rb");
        int len = 9;
	if (fp) {
		unsigned char *bytes;
		bytes = (unsigned char *) malloc (len * sizeof(unsigned char));
              
                if (fread(bytes, sizeof(unsigned char), len, fp)) {
                        mpz_import(tmp, len, 1, sizeof(unsigned char), 0, 0, bytes);
                }
                
                else {
                        printf("file read error\n");
                }
                
		fclose(fp);
		free(bytes);
	}
        
        else {
                printf("random number generation error\n");
        }
}

void guassian_poly(fmpz_poly_t poly, long d)
{
        if ( d == 0 )
		return;
	double tdvn = bgv_get_dvn();
	long a = (long)ceil(-10*tdvn);
	long b = (long)floor(+10*tdvn);
	long x, i;
	double p;
	mpz_t randseed;
	mpz_init(randseed);
	hcrypt_random(randseed);
	unsigned long int useed = mpz_get_ui(randseed);
	srand(useed);
	for( i = 0 ; i < d ; i++) {
                x = rand()%(b - a) + a;
                fmpz_poly_set_coeff_si(poly, i, x);
	}
	mpz_clear(randseed);
}

/*void unif_poly(fmpz_poly_t poly, fmpz_t space, long d)
{
        long i;
	mpz_t randseed;
	mpz_init(randseed);
	hcrypt_random(randseed);
	mpz_t rndnum, rndbd;
	fmpz_t rndfmpz;
        fmpz_init(rndfmpz);
	gmp_randstate_t gmpstate;
        
	mpz_init(rndnum);
	mpz_init(rndbd);
	fmpz_get_mpz(rndbd, space);
	
	gmp_randinit_default(gmpstate);
	gmp_randseed(gmpstate, randseed);
        
	for( i = 0 ; i < d ; i++ ) {
		mpz_urandomm(rndnum, gmpstate, rndbd);
		fmpz_set_mpz(rndfmpz, rndnum);
		fmpz_poly_set_coeff_fmpz(poly, i, rndfmpz);
	}
	mpz_clear(randseed);
	fmpz_clear(rndfmpz);
	gmp_randclear(gmpstate);
	mpz_clear(rndnum);
	mpz_clear(rndbd);
}*/

void unif_poly(fmpz_poly_t poly, fmpz_t space, long d)
{
        long i, x;
	mpz_t randseed;
	mpz_init(randseed);
	hcrypt_random(randseed);
	unsigned long int useed = mpz_get_ui(randseed);
	srand(useed);
	for( i = 0 ; i < d ; i++) {
                x = rand();
                fmpz_poly_set_coeff_si(poly, i, x);
	}
        fmpz_poly_scalar_smod_fmpz(poly, poly, space);
	mpz_clear(randseed);
}

