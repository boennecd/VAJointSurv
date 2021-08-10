#pragma once

#include "AADNumWrapper.h"
#include <numeric>

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
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of, bool trans){
        const size_t n = static_cast<size_t>(std::distance(af, al)), 
                     m = static_cast<size_t>(std::distance(xf, xl)) / n;
                     
        for(auto v = of; v != of + m; ++v)
            *v = 0;
        
        if(trans)
            for(size_t i = 0; i < m; ++i)
                for(size_t j = 0; j < n; ++j, ++xf)
                    of[i] += *xf * af[j];
        else
            for(size_t j = 0; j < n; ++j, ++af)
                for(size_t i = 0; i < m; ++i, ++xf)
                    of[i] += *xf * *af;
    }
};

template<class I1, class I2, class V2>
struct VecMatProdOp<I1, I2, Number, V2> {
    /// special case where the matrix is a Number type
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of, bool trans){
        Number::mat_vec_prod_TMat(xf, xl, af, al, of, trans);
    }
};

template<class I1, class I2, class V1>
struct VecMatProdOp<I1, I2, V1, Number> {
    /// special case where the vector is a Number type
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of, bool trans){
        Number::mat_vec_prod_TVec(xf, xl, af, al, of, trans);
    }
};

template<class I1, class I2>
struct VecMatProdOp<I1, I2, Number, Number> {
    /// special case where both iterators are for Numbers
    template<class I3>
    static void mat_vec_prod(I1 xf, I1 xl, I2 af, I2 al, I3 of, bool trans){
        Number::mat_vec_prod_identical(xf, xl, af, al, of, trans);
    }
};
} // namespace implementation

template<class I1, class I2, class I3>
void matVecProd(I1 xf, I1 xl, I2 af, I2 al, I3 of, bool trans){
    implementation::VecMatProdOp
    <I1, I2, it_value_type<I1>, it_value_type<I2> >::mat_vec_prod
    (xf, xl, af, al, of, trans);
}


// triangular matrix-vector product
namespace implementation {
template<class I1, class I2, class V1, class V2>
struct VecTriMatProdOp {
    /// the general case
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of, bool trans){
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
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of, bool trans){
        Number::trimat_vec_prod_Tmat(xf, af, al, of, trans);
    }
};

template<class I1, class I2, class V1>
struct VecTriMatProdOp<I1, I2, V1, Number>{
    /// special case where the vector is a Number type
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of, bool trans){
        Number::trimat_vec_prod_Tvec(xf, af, al, of, trans);
    }
};

template<class I1, class I2>
struct VecTriMatProdOp<I1, I2, Number, Number>{
    /// special case where both iterators are for Numbers
    template<class I3>
    static void trimat_vec_prod(I1 xf, I2 af, I2 al, I3 of, bool trans){
        Number::trimat_vec_prod_identical(xf, af, al, of, trans);
    }
};
} // namespace implementation

template<class I1, class I2, class I3>
void triMatVecProd(I1 xf, I2 af, I2 al, I3 of, bool trans){
    implementation::VecTriMatProdOp
    <I1, I2, it_value_type<I1>, it_value_type<I2> >::trimat_vec_prod
    (xf, af, al, of, trans);
}


} // namespace cfadd
