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

int main(int argc, char *args[])  /*bb, lam, lev*/
{
        fmpz_t b, lambda, level;
        fmpz_init(b);
        int flag = fmpz_set_str(b, args[1], 10);
        if(flag == -1){
                printf("invalid number b\n");
                exit(0);
        }
        long bb = fmpz_get_si(b);
        if(bb != 1 && bb != 0) {
                printf("the value of first parameter should be 0 or 1\n");
                exit(0);
        }
        
        fmpz_init(lambda);
        flag = fmpz_set_str(lambda, args[2], 10);
        if(flag == -1){
                printf("invalid number lambda\n");
                exit(0);
        }
        
        fmpz_init(level);
        flag = fmpz_set_str(level, args[3], 10);
        if(flag == -1){
                printf("invalid number level\n");
                exit(0);
        }
        
        long lam = fmpz_get_si(lambda);
        long lev = fmpz_get_si(level);
        fmpz_t mult;
        fmpz_init(mult);
        fmpz_mul(mult, lambda, level);
        long miu;
        miu = fmpz_flog_ui(mult, 2);
        param_node_t *param;
        param = param_node_init(param);
        long high;
        high = (lev + 1) * miu;
        gen_q(param->q, high);
        fmpz_t tmp;
        fmpz_init(tmp);
        fmpz_fdiv_q_si(tmp, param->q, bgv_get_bigb());
        long prod, d, mi;
        prod = lam * fmpz_flog_ui(tmp, 2);
        fmpz_set_si(tmp, prod);
        mi = fmpz_flog_ui(tmp, 2);
        prod = pow(2, mi);
        
        if(bb == 1) {
		param->n = prod;
		d = 1;
	}
	else {
		param->n = 1;
		d = prod;
	}
        param->bign = ceil((2 * param->n + 1) * fmpz_flog_ui(param->q, 2));
	param_node_t *r, *pn;
        
        r = param;
        long j;
        for(j = lev - 1 ; j >= 0 ; j--) {
                pn = e_setup((j+1)*miu, lam, bb, pn);
                r->next = pn;
                r = pn;
        }
        r->next = NULL;
        
        printf("%ld\n%ld\n", lev, d);
        r = param;
        while(r!=NULL) {
                fmpz_print(r->q);
                printf("\n%ld\n%ld\n", r->n, r->bign);
                r = r->next;
        }
        
        fmpz_clear(tmp);
        return 0;
}

/* output d {q, n, bign}
 */

