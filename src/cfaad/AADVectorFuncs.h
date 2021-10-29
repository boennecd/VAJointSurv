#pragma once

#include "AADNumWrapper.h"
#include <numeric>
#include <algorithm>

namespace cfaad {
// sum function
namespace implementation {
template<class I, class V>
struct VecSumOp {
    /// general case
    static V sum(I begin, I end){
        return std::accumulate(begin, end, V{0.});
    }
};

template<class I>
struct VecSumOp<I, Number> {
    /// sum over range T iterator
    static Number sum(I begin, I end){
        return Number::sum(begin, end);
    }
};
} // namespace implementation

template<class I>
it_value_type<I> sum(I begin, I end){
    return implementation::VecSumOp<I, it_value_type<I> >::sum(begin, end);
}

// dot product
namespace implementation {
template<class I1, class I2, class V1, class V2>
struct VecDotProdOp {
    using returnT = V1;
    /// general case
    static V1 dot_prodcut
    (I1 first1, I1 last1, I2 first2){
        return std::inner_product(first1, last1, first2, V1{0.});
    }
};

template<class I1, class I2, class V2>
struct VecDotProdOp<I1, I2, Number, V2> {
    using returnT = Number;
    /// dot product with one T and none T iterator
    static Number dot_prodcut
    (I1 first1, I1 last1, I2 first2){
        return Number::dot_product(first1, last1, first2);
    }
};

template<class I1, class I2, class V1>
struct VecDotProdOp<I1, I2, V1, Number> {
    using returnT = Number;
    /// dot product with one T and none T iterator
    static Number dot_prodcut
    (I1 first1, I1 last1, I2 first2){
        auto const n = std::distance(first1, last1);
        return Number::dot_product(first2, first2 + n, first1);
    }
};

template<class I1, class I2>
struct VecDotProdOp<I1, I2, Number, Number> {
    using returnT = Number;
    /// dot product with two T iterators
    static Number dot_prodcut
        (I1 first1, I1 last1, I2 first2){
        return Number::dot_product_identical(first1, last1, first2);
    }
};

template<class I1>
struct VecDotProdOp<I1, I1, Number, Number> {
    using returnT = Number;
    /// dot product with two T iterators of the same type
    static Number dot_prodcut
        (I1 first1, I1 last1, I1 first2){
        return Number::dot_product_identical_it(first1, last1, first2);
    }
};
} // namespace implementation

template<class I1, class I2>
typename implementation::VecDotProdOp
<I1, I2, it_value_type<I1>, it_value_type<I2> >::returnT
dotProd(I1 first1, I1 last1, I2 first2){
  return implementation::VecDotProdOp
      <I1, I2, it_value_type<I1>, it_value_type<I2> >::dot_prodcut
      (first1, last1, first2);
}

// matrix vector products
namespace implementation {
template<class I1, class I2, class V1, class V2>
struct VecMatProdOp {
    /// the general case
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                             const bool trans, const size_t other){
        if(trans){
            const size_t m{static_cast<size_t>(std::distance(af, al))},
                         n{other},
                       ldx{static_cast<size_t>(std::distance(xf, xl)) / n};

            for(size_t i = 0; i < n; ++i, xf += ldx - m){
                of[i] = 0;
                for(size_t j = 0; j < m; ++j, ++xf)
                    of[i] += *xf * af[j];
            }
            return;
        }

        const size_t m{other},
                     n{static_cast<size_t>(std::distance(af, al))},
                   ldx{static_cast<size_t>(std::distance(xf, xl)) / n};

        std::fill(of, of + m, 0);

        for(size_t j = 0; j < n; ++j, ++af, xf += ldx - m)
            for(size_t i = 0; i < m; ++i, ++xf)
                of[i] += *xf * *af;
    }
};

template<class I1, class I2, class V2>
struct VecMatProdOp<I1, I2, Number, V2> {
    /// special case where the matrix is a Number type
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                             const bool trans, const size_t other){
        Number::mat_vec_prod_TMat(xf, xl, af, al, of, trans, other);
    }
};

template<class I1, class I2, class V1>
struct VecMatProdOp<I1, I2, V1, Number> {
    /// special case where the vector is a Number type
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                             const bool trans, const size_t other){
        Number::mat_vec_prod_TVec(xf, xl, af, al, of, trans, other);
    }
};

template<class I1, class I2>
struct VecMatProdOp<I1, I2, Number, Number> {
    /// special case where both iterators are for Numbers
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                             const bool trans, const size_t other){
        Number::mat_vec_prod_identical(xf, xl, af, al, of, trans, other);
    }
};
} // namespace implementation

/// version where the leading dimension is the same as the matrix in the product
template<class I1, class I2, class I3>
void matVecProd(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                const bool trans){
    const size_t la{static_cast<size_t>(std::distance(af, al))},
              other{static_cast<size_t>(std::distance(xf, xl)) / la};

    implementation::VecMatProdOp
    <I1, I2, it_value_type<I1>, it_value_type<I2> >::mat_vec_prod
    (xf, xl, af, al, of, trans, other);
}

template<class I1, class I2, class I3>
void matVecProd(I1 xf, I1 xl, I2 af, I2 al, I3 of,
                const bool trans, const size_t other){
    implementation::VecMatProdOp
    <I1, I2, it_value_type<I1>, it_value_type<I2> >::mat_vec_prod
    (xf, xl, af, al, of, trans, other);
}


// triangular matrix-vector product
namespace implementation {
template<class I1, class I2, class V1, class V2>
struct VecTriMatProdOp {
    /// the general case
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of,
                                const bool trans){
        const size_t n = static_cast<size_t>(std::distance(af, al));

        for(auto v = of; v != of + n; ++v)
            *v = 0;

        if(trans)
            for(size_t j = 0; j < n; ++j)
                for(size_t i = 0; i <= j; ++i, ++xf)
                    of[j] += *xf * af[i];
        else
            for(size_t j = 0; j < n; ++j)
                for(size_t i = 0; i <= j; ++i, ++xf)
                    of[i] += *xf * af[j];
    }
};

template<class I1, class I2, class V2>
struct VecTriMatProdOp<I1, I2, Number, V2>{
    /// special case where the matrix is a Number type
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of,
                                const bool trans){
        Number::trimat_vec_prod_Tmat(xf, af, al, of, trans);
    }
};

template<class I1, class I2, class V1>
struct VecTriMatProdOp<I1, I2, V1, Number>{
    /// special case where the vector is a Number type
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of,
                                const bool trans){
        Number::trimat_vec_prod_Tvec(xf, af, al, of, trans);
    }
};

template<class I1, class I2>
struct VecTriMatProdOp<I1, I2, Number, Number>{
    /// special case where both iterators are for Numbers
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of,
                                const bool trans){
        Number::trimat_vec_prod_identical(xf, af, al, of, trans);
    }
};
} // namespace implementation

template<class I1, class I2, class I3>
void triMatVecProd(I1 xf, I2 af, I2 al, I3 of,
                   const bool trans){
    implementation::VecTriMatProdOp
    <I1, I2, it_value_type<I1>, it_value_type<I2> >::trimat_vec_prod
    (xf, af, al, of, trans);
}

#if AADLAPACK
// the quadratic form a^TB^{-1}a

namespace implementation {
template<class I1, class I2, class V1, class V2>
struct quadFormInvOp {
    /// the general case
    using returnT = double;
    static double quad_form_inv(I1 a, I2 B, const CholFactorization &chol){
        const size_t n{static_cast<size_t>(chol.n)};
        double *wk_mem{Number::tape->getWKMem(n)};

        // compute the result
        for(size_t i = 0; i < n; ++i, ++a)
            wk_mem[i] = *a;
        chol.solveU(wk_mem, true); // U^{-T}a
        double res{0};
        for(size_t i = 0; i < n; ++i, ++wk_mem)
            res += *wk_mem * *wk_mem;

        return res;
    }
};

template<class I1, class I2, class V2>
struct quadFormInvOp<I1, I2, Number, V2> {
    /// special case where the vector is a Number iterator
    using returnT = Number;
    static Number quad_form_inv(I1 a, I2 B, const CholFactorization &chol){
        return Number::quadFormInv_TVec(a, B, chol);
    }
};

template<class I1, class I2, class V1>
struct quadFormInvOp<I1, I2, V1, Number> {
    /// special case where the matrix is a Number iterator
    using returnT = Number;
    static Number quad_form_inv(I1 a, I2 B, const CholFactorization &chol){
        return Number::quadFormInv_TMat(a, B, chol);
    }
};

template<class I1, class I2>
struct quadFormInvOp<I1, I2, Number, Number> {
    /// special case where the both iterators are for Numbers
    using returnT = Number;
    static Number quad_form_inv(I1 a, I2 B, const CholFactorization &chol){
        return Number::quadFormInv_identical(a, B, chol);
    }
};
} // namespace implementation

template<class I1, class I2>
typename implementation::quadFormInvOp
<I1, I2, it_value_type<I1>, it_value_type<I2> >::returnT
quadFormInv(I1 a, I2 B, const CholFactorization &chol){
    return implementation::quadFormInvOp
        <I1, I2, it_value_type<I1>, it_value_type<I2> >
        ::quad_form_inv(a, B, chol);
}

// the trace tr(A^{-1}B) for two quadratic matrices

namespace implementation {
template<class I1, class I2, class V1, class V2>
struct trInvMatMatOp {
    /// the general case
    using returnT = double;
    static double tr_invmat_mat(I1 A, I2 B, const CholFactorization &chol){
        const size_t n{static_cast<size_t>(chol.n)},
                    nn = n * n;
        double *wk_mem{Number::tape->getWKMem(nn)};

        // compute the result
        for(size_t i = 0; i < nn; ++i, ++B)
            wk_mem[i] = *B;
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem + i * n); // A^{-1}B

        double res{0};
        for(size_t i = 0; i < n; ++i, wk_mem += n + 1)
            res += *wk_mem;

        return res;
    }
};

template<class I1, class I2, class V2>
struct trInvMatMatOp<I1, I2, Number, V2> {
    using returnT = Number;
    /// special case where the first iterator is to Ts
    static Number tr_invmat_mat(I1 A, I2 B, const CholFactorization &chol){
        return Number::trInvMatMat_first(A, B, chol);
    }
};

template<class I1, class I2, class V1>
struct trInvMatMatOp<I1, I2, V1, Number> {
    using returnT = Number;
    /// special case where the second iterator is to Ts
    static Number tr_invmat_mat(I1 A, I2 B, const CholFactorization &chol){
        return Number::trInvMatMat_second(A, B, chol);
    }
};

template<class I1, class I2>
struct trInvMatMatOp<I1, I2, Number, Number> {
    using returnT = Number;
    /// special case where the both iterators are to Ts
    static Number tr_invmat_mat(I1 A, I2 B, const CholFactorization &chol){
        return Number::trInvMatMat_identical(A, B, chol);
    }
};

} // namespace implementation

template<class I1, class I2>
typename implementation::trInvMatMatOp
<I1, I2, it_value_type<I1>, it_value_type<I2> >::returnT
trInvMatMat(I1 A, I2 B, const CholFactorization &chol){
    return implementation::trInvMatMatOp
        <I1, I2, it_value_type<I1>, it_value_type<I2> >
        ::tr_invmat_mat(A, B, chol);
}

// the log determinant

namespace implementation {
template<class I, class V>
struct logDeterOp {
    using returnT = double;
    /// the general case
    static double log_deter(I begin, const CholFactorization &chol){
        return log(chol.determinant());
    }
};

template<class I>
struct logDeterOp<I, Number> {
    using returnT = Number;
    /// special case where the iterator is for Numbers
    static Number log_deter(I begin, const CholFactorization &chol){
        return Number::log_deter(begin, chol);
    }
};
} // namespace implementation

template<class I>
typename implementation::logDeterOp<I, it_value_type<I> >::returnT
logDeter(I begin, const CholFactorization &chol){
    return implementation::logDeterOp<I, it_value_type<I> >
        ::log_deter(begin, chol);
}

// quadratic form

namespace implementation {
template<class I1, class I2, class V1, class V2>
struct quadFormSymOP {
    /// the general case
    using returnT = double;
    static double quad_form_sym(I1 x, I2 yf, I2 ye){
        const size_t n{static_cast<size_t>(std::distance(yf, ye))};

        // compute the result
        double res{};
        for(size_t j = 0; j < n; ++j, x += n){
            for(size_t i = 0; i < j; ++i)
                res += x[i] * yf[i] * yf[j];
            res += .5 * x[j] * yf[j] * yf[j];
        }

        return 2 * res;
    }
};

template<class I1, class I2, class V2>
struct quadFormSymOP<I1, I2, Number, V2> {
    using returnT = Number;
    static Number quad_form_sym(I1 x, I2 yf, I2 ye){
        return Number::quad_form_sym_mat(x, yf, ye);
    }
};

template<class I1, class I2, class V1>
struct quadFormSymOP<I1, I2, V1, Number> {
    using returnT = Number;
    static Number quad_form_sym(I1 x, I2 yf, I2 ye){
        return Number::quad_form_sym_vec(x, yf, ye);
    }
};

template<class I1, class I2>
struct quadFormSymOP<I1, I2, Number, Number> {
    using returnT = Number;
    static Number quad_form_sym(I1 x, I2 yf, I2 ye){
        return Number::quad_form_sym_both(x, yf, ye);
    }
};

} // namespace implementation

template<class I1, class I2>
typename implementation::quadFormSymOP
<I1, I2, it_value_type<I1>, it_value_type<I2> >::returnT
quadFormSym(I1 x, I2 yf, I2 ye){
    return implementation::quadFormSymOP
      <I1, I2, it_value_type<I1>, it_value_type<I2> >
        ::quad_form_sym(x, yf, ye);
}

#endif // if AADLAPACK

} // namespace cfaad
