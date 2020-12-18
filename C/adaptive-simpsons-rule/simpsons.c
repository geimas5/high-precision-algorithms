#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

const int bits = 256;


mpfr_t from;
mpfr_t to;
mpfr_t result;
mpfr_t eps;
mpfr_t eps_x_15;
mpfr_t eps_div_2;

//tmp for use within approximations.
mpfr_t f_tmp_1;

void _quad_simpsons_mem(mpfr_ptr r_m, mpfr_ptr r_fm, mpfr_ptr r_simp, void (*f)(mpfr_ptr, mpfr_ptr), mpfr_ptr a, mpfr_ptr fa, mpfr_ptr b, mpfr_ptr fb) {
    mpfr_t tmp1;
    mpfr_init2 (tmp1, bits);

    mpfr_add(r_m, a, b, MPFR_RNDD);
    mpfr_div_si(r_m, r_m, 2, MPFR_RNDD);

    f(r_fm, r_m);

    mpfr_mul_si(tmp1, r_fm, 4, MPFR_RNDD);
    mpfr_add( tmp1, tmp1, fb, MPFR_RNDD );
    mpfr_add( tmp1, tmp1, fa, MPFR_RNDD );

    mpfr_sub(r_simp, b, a, MPFR_RNDD);
    mpfr_div_si(r_simp, r_simp, 6, MPFR_RNDD);

    mpfr_mul( r_simp, tmp1, r_simp, MPFR_RNDD );

    mpfr_clear(tmp1);
}

void _quad_asr(mpfr_ptr r, void (*f)(mpfr_ptr, mpfr_ptr), mpfr_ptr a, mpfr_ptr fa, mpfr_ptr b, mpfr_ptr fb, mpfr_ptr eps, mpfr_ptr whole, mpfr_ptr m, mpfr_ptr fm) {

    mpfr_t lt_m, lt_fm, lt_simp, rt_m, rt_fm, rt_simp, delta, tmp;
    mpfr_init2 (lt_m, bits);
    mpfr_init2 (lt_fm, bits);
    mpfr_init2 (lt_simp, bits);
    mpfr_init2 (rt_m, bits);
    mpfr_init2 (rt_fm, bits);
    mpfr_init2 (rt_simp, bits);
    mpfr_init2 (delta, bits);
    mpfr_init2 (tmp, bits);

    _quad_simpsons_mem(lt_m, lt_fm, lt_simp, f, a, fa, m, fm);
    _quad_simpsons_mem(rt_m, rt_fm, rt_simp, f, m, fm, b, fb);

    mpfr_add( delta, lt_simp, rt_simp, MPFR_RNDD );
    mpfr_sub( delta, delta, whole, MPFR_RNDD );
    mpfr_abs( delta, delta, MPFR_RNDD );

    if (mpfr_lessequal_p(delta, eps_x_15)) {
        mpfr_div_si( r, delta, 15, MPFR_RNDD );
        mpfr_add( r, r, lt_simp, MPFR_RNDD );
        mpfr_add( r, r, rt_simp, MPFR_RNDD );
    }
    else{
        _quad_asr(r, f, a, fa, m, fm, eps_div_2, lt_simp, lt_m, lt_fm);
        _quad_asr(tmp, f, m, fm, b, fb, eps_div_2, rt_simp, rt_m, rt_fm);
        mpfr_add(r, r, tmp, MPFR_RNDD );
    }

    mpfr_clear (lt_m);
    mpfr_clear (lt_fm);
    mpfr_clear (lt_simp);
    mpfr_clear (rt_m);
    mpfr_clear (rt_fm);
    mpfr_clear (rt_simp);
    mpfr_clear (delta);
    mpfr_clear (tmp);
}

void quad_asr(mpfr_ptr r, void (*f)(mpfr_ptr, mpfr_ptr), mpfr_ptr a, mpfr_ptr b, mpfr_ptr eps) {
    mpfr_t fa, fb, m, fm, simp;
    mpfr_init2 (fa, bits);
    mpfr_init2 (fb, bits);
    mpfr_init2 (m, bits);
    mpfr_init2 (fm, bits);
    mpfr_init2 (simp, bits);

    f(fa, a);
    f(fb, b);

    _quad_simpsons_mem(m, fm, simp, f, a, fa, b, fb);

    _quad_asr(r, f, a, fa, b, fb, eps, simp, m, fm);

    mpfr_clear (m);
    mpfr_clear (fm);
    mpfr_clear (simp);
}
 
void sinx(mpfr_ptr rop, mpfr_ptr x ){
    mpfr_sin( rop, x, MPFR_RNDD );
}



int main (void)
{
    //integrate_test();
    //return 0;

  mpfr_init2 (result, bits);
  mpfr_init2 (from, bits);
  mpfr_init2 (to, bits);
  mpfr_init2 (eps, bits);
  mpfr_init2 (eps_x_15, bits);
  mpfr_init2 (eps_div_2, bits);

  mpfr_set_str (from, "3", 10, MPFR_RNDD);
  mpfr_set_str (to, "4", 10, MPFR_RNDD);
  mpfr_set_d (eps, 1e-10, MPFR_RNDD);


  mpfr_mul_si( eps_x_15, eps, 15, MPFR_RNDD );
  mpfr_div_si( eps_div_2, eps, 2, MPFR_RNDD );

  quad_asr( result, sinx, from, to, eps );
  mpfr_printf( "Sin from %.RF to %.0RF = %.10RF\n", from, to, result );

  mpfr_clear(result);
  mpfr_clear(from);
  mpfr_clear(to);
  mpfr_clear(eps);
  mpfr_clear(eps_x_15);
  mpfr_clear(eps_div_2);

  mpfr_free_cache ();
  return 0;
}