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

char str[100000];

void print(fmpz_poly_mat_t poly);

int main(int argc, char *args[])/*setup.txt,sk.txt{row, col, poly}*/
{
	FILE *fp;
        
        if((fp = fopen(args[1], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
                
        fmpz_t tmp;
        fmpz_init(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        
        long lev, d, i, j;
        lev = fmpz_get_si(tmp);
        fgets(str, 100, fp);
        fmpz_set_str(tmp, str, 10);
        d = fmpz_get_si(tmp);
        
        fmpz_poly_t fx;
        fmpz_poly_init(fx);
        fmpz_poly_set_coeff_si(fx, 0, 1);
        fmpz_poly_set_coeff_si(fx, d, 1);
        
        param_node_t *ph, *pr, *ps, *param, *pam;
        ph = param_node_init(ph);

        ps = ph;
        for( i = 0 ; i <= lev; i++ ) {
                pr = param_node_init(pr);
                fgets(str, 100, fp);
                fmpz_set_str(pr->q, str, 10);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->n = fmpz_get_si(tmp);
                fgets(str, 100, fp);
                fmpz_set_str(tmp, str, 10);
                pr->bign = fmpz_get_si(tmp);
                ps->next = pr;
                ps = pr;
        }
        ps->next = NULL;
        
        fclose(fp);
        long row, col;

        if((fp = fopen(args[2], "r")) == NULL)
        {
                printf("file read error\n");
                exit(0);
        }
        
        long l;
        sk_node_t *sh, *ss, *sr;
        sh = (sk_node_t *)malloc(sizeof(sk_node_t));
        
        ss = sh;
        for(l = 0; l <= lev ; l++){
                sr = (sk_node_t *)malloc(sizeof(sk_node_t));
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                row = fmpz_get_si(tmp);                
                
                fgets(str, 30, fp);
                fmpz_set_str(tmp, str,10);
                col = fmpz_get_si(tmp);
                
                fmpz_poly_mat_init(sr->sk, row, col);
                
                for( i = 0 ; i < row ; i++) {
                        for(j = 0; j < col ; j++) {
                                fgets(str, 100000, fp);
                                fmpz_poly_set_str(fmpz_poly_mat_entry(sr->sk, i, j), str);
                        }
                }
                ss->next = sr;
                ss = sr;
        }
        ss->next = NULL;
        
        fclose(fp);
        
        param = ph->next;
        pam = param->next;
        
        sh = sh->next;
        fmpz_poly_mat_t pka,pkb;
        fmpz_poly_mat_init(pkb, 1, 1);
        fmpz_poly_mat_init(pka, param->bign, 1 + (param->n));
        e_pkeygen(pka, param, sh->sk, d, fx);
        print(pka);
        print(pkb);
        fmpz_poly_mat_clear(pka);
        fmpz_poly_mat_clear(pkb);
        ss = sh;
        sr = ss->next;
        
        fmpz_poly_mat_t s1, tensor;
	long row1, row2, len, llog;
        
        for( i = lev ; i > 0 ; i-- ) {
                llog = fmpz_clog_ui(pam->q, 2);
		fmpz_poly_mat_init(pka, pam->bign, 1 + (pam->n));

                e_pkeygen(pka, pam, sr->sk, d, fx);
                
                row1 = fmpz_poly_mat_nrows(ss->sk);
        	row2 = row1 * row1;
        	fmpz_poly_mat_init(tensor, row2, 1);
                vec_tensor(tensor, ss->sk, param->q, fx);
                
                len = fmpz_clog_ui(param->q, 2);
		row2 = row2 * len;
                fmpz_poly_mat_init(s1, row2, 1);
                bitdecomp(s1, tensor, param->q, d);
                
                row1 = fmpz_poly_mat_nrows(s1) * llog;
		row2 = fmpz_poly_mat_nrows(sr->sk);
		fmpz_poly_mat_init(pkb, row1, row2);
                
                switchkeygen(pkb, s1, sr->sk, pam->q, d, fx);
                
                fmpz_poly_mat_clear(s1);
                fmpz_poly_mat_clear(tensor);
                pam = pam->next;
                param = param->next;
                ss = sr;
                sr = sr->next;
                print(pka);
                print(pkb);
                fmpz_poly_mat_clear(pka);
                fmpz_poly_mat_clear(pkb);
        }
        
        return 0;
}

void print(fmpz_poly_mat_t poly)
{
        long row, col, i, j;
        row = fmpz_poly_mat_nrows(poly);
        col = fmpz_poly_mat_ncols(poly);
        printf("%ld\n%ld\n", row, col);
        for( i = 0 ; i < row ; i++ ) {
                for( j = 0 ; j < col ; j++ ) {
                        fmpz_poly_print(fmpz_poly_mat_entry(poly, i, j));
                        printf("\n");
                }
        }
}