#pragma once

#include <iterator>
#include <type_traits>
#include "AADLp.h"
#include <algorithm>

namespace cfaad {

template<class T>
struct vectorOps {
    /// computes the sum
    template<class I>
    static T sum(I begin, I end){
        static_assert(is_it_value_type<I, T>::value, "Iterator is not to Ts");

        T res;
        res.createNode(static_cast<size_t>(std::distance(begin, end)));
        res.myValue = 0;
        for(size_t i = 0; begin != end; ++begin, ++i){
            res.myValue += begin->value();
            res.setpDerivatives(i, 1.);
            res.setpAdjPtrs(i, *begin);
        }

        return res;
    }

    /// computes the dot product with T and none T iterator
    template<class I1, class I2>
    static T dot_product(I1 f1, I1 l1, I2 f2){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "Second iterator is to Ts");

        T res;
        res.createNode(static_cast<size_t>(std::distance(f1, l1)));
        res.myValue = 0;
        for(size_t i = 0; f1 != l1; ++f1, ++f2, ++i){
          res.myValue += f1->value() * *f2;
          res.setpDerivatives(i, *f2);
          res.setpAdjPtrs(i, *f1);
        }

        return res;
    }

    /// computes the dot product with T iterators
    template<class I1, class I2>
    static T dot_product_identical(I1 f1, I1 l1, I2 f2){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
                      
        T res;
        const size_t n{static_cast<size_t>(std::distance(f1, l1))};
        res.createNode(2 * n);
        res.myValue = 0;
        for(size_t i = 0; f1 != l1; ++f1, ++f2, ++i){
            res.myValue += f1->value() * f2->value();
            res.setpDerivatives(i, f2->value());
            res.setpAdjPtrs(i, *f1);
            res.setpDerivatives(i + n, f1->value());
            res.setpAdjPtrs(i + n, *f2);
        }

        return res;
    }
    
    /// special case of with two T iterators of the same type
    template<class I1>
    static T dot_product_identical_it(I1 f1, I1 l1, I1 f2){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
                      
        if(f1 != f2)
            return dot_product_identical(f1, l1, f2);
                      
        // special case where the iterators are the same
        T res;
        const size_t n{static_cast<size_t>(std::distance(f1, l1))};        
        res.createNode(n);
        res.myValue = 0;
        for(size_t i = 0; f1 != l1; ++f1, ++i){
            res.myValue += f1->value() * f1->value();
            res.setpDerivatives(i, 2 * f1->value());
            res.setpAdjPtrs(i, *f1);
        }

        return res;
    }
    
    /*
     * computes the matrix vector product X.a where X is a m x n matrix and 
     * a is a n vector. The result is stored in the last argument which needs
     * to be iterator with space for at least m elements. The matrix is in 
     * column-major order.
     * 
     * This is the version where the matrix is a T type whereas the vector is 
     * for non-Ts. The trans argument sets whether it is the X^T rather than X.
     * The trans argument determines m when trans is FALSE or n when trans is 
     * TRUE.
     */
    template<class I1, class I2, class I3>
    static void mat_vec_prod_TMat
    (I1 xf, I1 xl, I2 af, I2 al, I3 of, const bool trans, const size_t other){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "Second iterator is to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                     
        if(trans){
            const size_t m{static_cast<size_t>(std::distance(af, al))}, 
                         n{other},
                       ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
            
            for(size_t i = 0; i < n; ++i, ++of, xf += ldx - m){
                of->createNode(m);
                of->myValue = 0;
                for(size_t j = 0; j < m; ++j, ++xf){
                    of->myValue += xf->value() * af[j];
                    of->setpDerivatives(j, af[j]);
                    of->setpAdjPtrs(j, *xf);
                }
            }
            
            return;
        }
        
        const size_t m{other},
                     n{static_cast<size_t>(std::distance(af, al))},
                   ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
                     
        for(auto v = of; v != of + m; ++v){
            v->createNode(n);
            v->myValue = 0;
        }
        
        for(size_t j = 0; j < n; ++j, ++af, xf += ldx - m)
            for(size_t i = 0; i < m; ++i, ++xf){
                of[i].myValue += xf->value() * *af;
                of[i].setpDerivatives(j, *af);
                of[i].setpAdjPtrs(j, *xf);
            }
    }
    
    /**
     * the same as mat_vec_prod_TMat but where the first argument is for 
     * non-Ts and the second argument is for Ts.
     */
    template<class I1, class I2, class I3>
    static void mat_vec_prod_TVec
    (I1 xf, I1 xl, I2 af, I2 al, I3 of, const bool trans, const size_t other){
        static_assert(!is_it_value_type<I1, T>::value,
                      "First iterator is to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                     
        if(trans){
            const size_t m{static_cast<size_t>(std::distance(af, al))}, 
                         n{other},
                       ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
            
            for(size_t i = 0; i < n; ++i, ++of, xf += ldx - m){
                of->createNode(m);
                of->myValue = 0;
                for(size_t j = 0; j < m; ++j, ++xf){
                    of->myValue += *xf * af[j].value();
                    of->setpDerivatives(j, *xf);
                    of->setpAdjPtrs(j, af[j]);
                }
            }
            
            return;
        }
        
        const size_t m{other},
                     n{static_cast<size_t>(std::distance(af, al))},
                   ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
                     
        for(auto v = of; v != of + m; ++v){
            v->createNode(n);
            v->myValue = 0;
        }
        
        for(size_t j = 0; j < n; ++j, ++af, xf += ldx - m)
            for(size_t i = 0; i < m; ++i, ++xf){
                of[i].myValue += *xf  * af->value();
                of[i].setpDerivatives(j, *xf);
                of[i].setpAdjPtrs(j, *af);
            }
    }
    
    /// the same as mat_vec_prod_TMat but where both iterators are for Ts.
    template<class I1, class I2, class I3>
    static void mat_vec_prod_identical
    (I1 xf, I1 xl, I2 af, I2 al, I3 of, const bool trans, const size_t other){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                     
        if(trans){
            const size_t m{static_cast<size_t>(std::distance(af, al))}, 
                         n{other},
                       ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
            
            for(size_t i = 0; i < n; ++i, ++of, xf += ldx - m){
                of->createNode(2 * m);
                of->myValue = 0;
                for(size_t j = 0; j < m; ++j, ++xf){
                    of->myValue += xf->value() * af[j].value();
                    of->setpDerivatives(j, xf->value());
                    of->setpAdjPtrs(j, af[j]);
                    of->setpDerivatives(j + m, af[j].value());
                    of->setpAdjPtrs(j + m, *xf);
                }
            }
            
            return;
        }
        
        const size_t m{other},
                     n{static_cast<size_t>(std::distance(af, al))},
                   ldx{static_cast<size_t>(std::distance(xf, xl)) / n};
                     
        for(auto v = of; v != of + m; ++v){
            v->createNode(2 * n);
            v->myValue = 0;
        }

        for(size_t j = 0; j < n; ++j, ++af, xf += ldx - m)
            for(size_t i = 0; i < m; ++i, ++xf){
                of[i].myValue += xf->value() * af->value();
                of[i].setpDerivatives(j, xf->value());
                of[i].setpAdjPtrs(j, *af);
                of[i].setpDerivatives(j + n, af->value());
                of[i].setpAdjPtrs(j + n, *xf);
            }
    }
    
    /*
     * computes the matrix vector product X.a where X is a n x n triangular 
     * matrix and a is a n vector. The result is stored in the last argument 
     * which needs to be iterator with space for at least n elements. The matrix 
     * is in column-major order and only the non-zero elements are passed.
     * 
     * This is the version where the matrix is a T type whereas the vector is 
     * for non-Ts. The last argument sets whether it is the X^T rather than X.
     */
    template<class I1, class I2, class I3>
    static void trimat_vec_prod_Tmat(I1 xf, I2 af, I2 al, I3 of, 
                                     const bool trans){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "Second iterator is to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                      
        const size_t n = static_cast<size_t>(std::distance(af, al));
        
        if(trans)
            for(size_t j = 0; j < n; ++j, ++of){
                of->createNode(j + 1);
                of->myValue = 0;
                for(size_t i = 0; i <= j; ++i, ++xf){
                    of->myValue += xf->value() * af[i];
                    of->setpDerivatives(i, af[i]);
                    of->setpAdjPtrs(i, *xf);
                }
            }
        else 
            for(size_t j = 0; j < n; ++j){
                of[j].createNode(n - j);
                of[j].myValue = 0;
                size_t p_idx{j};
                for(size_t i = 0; i <= j; ++i, ++xf, --p_idx){
                    of[i].myValue += xf->value() * af[j];
                    of[i].setpDerivatives(p_idx, af[j]);
                    of[i].setpAdjPtrs(p_idx, *xf);
                }
            }
    }
    
    /**
     * the same as trimat_vec_prod_Tmat but where the first argument is for 
     * non-Ts and the second argument is for Ts.
     */
    template<class I1, class I2, class I3>
    static void trimat_vec_prod_Tvec(I1 xf, I2 af, I2 al, I3 of, 
                                     const bool trans){
        static_assert(!is_it_value_type<I1, T>::value,
                      "First iterator is to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                      
        const size_t n = static_cast<size_t>(std::distance(af, al));
        
        if(trans)
            for(size_t j = 0; j < n; ++j, ++of){
                of->createNode(j + 1);
                of->myValue = 0;
                for(size_t i = 0; i <= j; ++i, ++xf){
                    of->myValue += *xf * af[i].value();
                    of->setpDerivatives(i, *xf);
                    of->setpAdjPtrs(i, af[i]);
                }
            }
        else 
            for(size_t j = 0; j < n; ++j){
                of[j].createNode(n - j);
                of[j].myValue = 0;
                size_t p_idx{j};
                for(size_t i = 0; i <= j; ++i, ++xf, --p_idx){
                    of[i].myValue += *xf * af[j].value();
                    of[i].setpDerivatives(p_idx, *xf);
                    of[i].setpAdjPtrs(p_idx, af[j]);
                }
            }
    }
    
    /// the same as trimat_vec_prod_Tmat but where both iterators are for Ts.
    template<class I1, class I2, class I3>
    static void trimat_vec_prod_identical(I1 xf, I2 af, I2 al, I3 of, 
                                          const bool trans){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        static_assert(is_it_value_type<I3, T>::value,
                      "Third iterators is not to Ts");
                      
        const size_t n = static_cast<size_t>(std::distance(af, al));
        
        if(trans)
            for(size_t j = 0; j < n; ++j, ++of){
                of->createNode(2 * (j + 1));
                of->myValue = 0;
                for(size_t i = 0; i <= j; ++i, ++xf){
                    of->myValue += xf->value() * af[i].value();
                    of->setpDerivatives(i, xf->value());
                    of->setpAdjPtrs(i, af[i]);
                    of->setpDerivatives(i + j + 1, af[i].value());
                    of->setpAdjPtrs(i + j + 1, *xf);
                }
            }
        else 
            for(size_t j = 0; j < n; ++j){
                of[j].createNode(2 * (n - j));
                of[j].myValue = 0;
                size_t p_idx{j};
                for(size_t i = 0; i <= j; ++i, ++xf, --p_idx){
                    of[i].myValue += xf->value() * af[j].value();
                    of[i].setpDerivatives(p_idx, xf->value());
                    of[i].setpAdjPtrs(p_idx, af[j]);
                    of[i].setpDerivatives(p_idx + n - i, af[j].value());
                    of[i].setpAdjPtrs(p_idx + n - i, *xf);
                }
            }
    }
    
    template<class I1, class I2>
    static T quad_form_sym_mat(I1 x, I2 yf, I2 ye){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "First iterator is to Ts");
                      
        const size_t n{static_cast<size_t>(std::distance(yf, ye))};
        T res;
        res.createNode(n * n);
        res.myValue = 0;
        
        for(size_t j = 0; j < n; ++j){
            for(size_t i = 0; i < j; ++i){
                double const y_prod{yf[i] * yf[j]};
                res.myValue += x[i + j * n].value() * y_prod;
                
                res.setpAdjPtrs    (i + j * n, x[i + j * n]);
                res.setpDerivatives(i + j * n, y_prod);
                res.setpAdjPtrs    (j + i * n, x[j + i * n]);
                res.setpDerivatives(j + i * n, y_prod);
                
            }
            
            double const y_prod{yf[j] * yf[j]};
            res.myValue += .5 * x[j + j * n].value() * y_prod;
            res.setpAdjPtrs    (j + j * n, x[j + j * n]);
            res.setpDerivatives(j + j * n, y_prod);
        }
        
        res.myValue *= 2;
        return res;
    }
    
    template<class I1, class I2>
    static T quad_form_sym_vec(I1 x, I2 yf, I2 ye){
        static_assert(!is_it_value_type<I1, T>::value,
                      "First iterator is to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "First iterator is not to Ts");
                      
        const size_t n{static_cast<size_t>(std::distance(yf, ye))};
        T res;
        res.createNode(n);
        res.myValue = 0;
        double *derivs{T::tape->getWKMem(n)};
        std::fill(derivs, derivs + n, 0);
        
        for(size_t j = 0; j < n; ++j){
            for(size_t i = 0; i < j; ++i){
                double const y_prod{yf[i].value() * yf[j].value()};
                res.myValue += x[i + j * n] * y_prod;
                derivs[i] += x[i + j * n] * yf[j].value();
                derivs[j] += x[i + j * n] * yf[i].value();
            }
            
            double const y_prod{yf[j].value() * yf[j].value()};
            res.myValue += .5 * x[j + j * n] * y_prod;
            derivs[j] += x[j + j * n] *  yf[j].value();
        }
        
        for(size_t i = 0; i < n; ++i){
            res.setpAdjPtrs    (i, yf[i]);
            res.setpDerivatives(i, 2 * derivs[i]);
        }
        
        res.myValue *= 2;
        return res;
    }
    
    template<class I1, class I2>
    static T quad_form_sym_both(I1 x, I2 yf, I2 ye){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "First iterator is not to Ts");
                      
        const size_t n{static_cast<size_t>(std::distance(yf, ye))};
        T res;
        res.createNode(n * (n + 1));
        res.myValue = 0;
        double *derivs{T::tape->getWKMem(n)};
        std::fill(derivs, derivs + n, 0);
        
        for(size_t j = 0; j < n; ++j){
            for(size_t i = 0; i < j; ++i){
                double const y_prod{yf[i].value() * yf[j].value()};
                res.myValue += x[i + j * n].value() * y_prod;
                
                res.setpAdjPtrs    (i + j * n, x[i + j * n]);
                res.setpDerivatives(i + j * n, y_prod);
                res.setpAdjPtrs    (j + i * n, x[j + i * n]);
                res.setpDerivatives(j + i * n, y_prod);
                
                derivs[i] += x[i + j * n].value() * yf[j].value();
                derivs[j] += x[i + j * n].value() * yf[i].value();
            }
            
            double const y_prod{yf[j].value() * yf[j].value()};
            res.myValue += .5 * x[j + j * n].value() * y_prod;
            res.setpAdjPtrs    (j + j * n, x[j + j * n]);
            res.setpDerivatives(j + j * n, y_prod);
            derivs[j] += x[j + j * n].value() *  yf[j].value();
        }
        
        for(size_t i = 0; i < n; ++i){
            res.setpAdjPtrs    (i + n * n, yf[i]);
            res.setpDerivatives(i + n * n, 2 * derivs[i]);
        }
        
        res.myValue *= 2;
        return res;
    }
    
#if AADLAPACK

    /**
     * computes the quadratic form a^TB^{-1}b where B is a symmetric positive 
     * definite matrix where a Cholesky factorization has already been computed.
     *
     * Here, B is an iterator for Ts while a is an iterator for doubles.
     */
    template<class I1, class I2>
    static T quadFormInv_TMat(I1 a, I2 B, const CholFactorization &chol){
        static_assert(!is_it_value_type<I1, T>::value,
                      "First iterator is to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");

        const size_t n = static_cast<size_t>(chol.n);
        double *wk_mem = T::tape->getWKMem(n);
        
        // compute the result
        using I1_diff_type = typename std::iterator_traits<I1>::difference_type;
        std::copy(a, a + I1_diff_type(n), wk_mem);
        chol.solveU(wk_mem, true); // U^{-T}a
        double out_value{};
        for(size_t i = 0; i < n; ++i)
            out_value += wk_mem[i] * wk_mem[i];
            
        /* compute the partial derivatives. These are given by 
         * -B^{-1}aa^\topB^{-1} */
        chol.solveU(wk_mem, false);
        T res;
        res.myValue = out_value;
        res.createNode(n * n);
        for(size_t j = 0; j < n; ++j)
            for(size_t i = 0; i < n; ++i, ++B){
                res.setpDerivatives(i + j * n, -wk_mem[i] * wk_mem[j]);
                res.setpAdjPtrs(i + j * n, *B);
            }
        
        return res;
    }
    
    /// same as quadFormInv_TMat but where the vector is to Ts
    template<class I1, class I2>
    static T quadFormInv_TVec(I1 a, I2 B, const CholFactorization &chol){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "Second iterator is to Ts");

        const size_t n = static_cast<size_t>(chol.n);
        double *wk_mem = T::tape->getWKMem(n);
        
        // compute the result
        {
            I1 ai = a;
            for(size_t i = 0; i < n; ++i, ++ai)
                wk_mem[i] = ai->value();
        }
        chol.solveU(wk_mem, true); // U^{-T}a
        double out_value{};
        for(size_t i = 0; i < n; ++i)
            out_value += wk_mem[i] * wk_mem[i];
            
        // compute the partial derivatives. These are given by 2 * B^{-1}a
        chol.solveU(wk_mem, false);
        T res;
        res.myValue = out_value;
        res.createNode(n);
        for(size_t i = 0; i < n; ++i, ++a){
            res.setpDerivatives(i, 2 * wk_mem[i]);
            res.setpAdjPtrs(i, *a);
        }
        
        return res;
    }
    
    /// same as quadFormInv_TMat but with two iterators to Ts
    template<class I1, class I2>
    static T quadFormInv_identical(I1 a, I2 B, const CholFactorization &chol){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");

        const size_t n = static_cast<size_t>(chol.n);
        double *wk_mem = T::tape->getWKMem(n);
        
        // compute the result
        {
            I1 ai = a;
            for(size_t i = 0; i < n; ++i, ++ai)
                wk_mem[i] = ai->value();
        }
        chol.solveU(wk_mem, true); // U^{-T}a
        double out_value{};
        for(size_t i = 0; i < n; ++i)
            out_value += wk_mem[i] * wk_mem[i];
            
        // compute the partial derivatives
        chol.solveU(wk_mem, false);
        T res;
        res.myValue = out_value;
        res.createNode(n * (n + 1));
        for(size_t j = 0; j < n; ++j, ++a){
            res.setpDerivatives(j, 2 * wk_mem[j]);
            res.setpAdjPtrs(j, *a);
            for(size_t i = 0; i < n; ++i, ++B){
                res.setpDerivatives(i + (j + 1) * n, -wk_mem[i] * wk_mem[j]);
                res.setpAdjPtrs(i + (j + 1) * n, *B);
            }
        }
        
        return res;
    }
    
    /**
     * Computes tr(A^{-1}B) for two symmetric square matrices using a 
     * pre-computed Cholesky decomposition. This is the version where the A 
     * iterator is for Ts and the B iterator is for non-Ts */
    template<class I1, class I2>
    static T trInvMatMat_first(I1 A, I2 B, const CholFactorization &chol){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(!is_it_value_type<I2, T>::value,
                      "Second iterator is to Ts");
        
        const size_t n{static_cast<size_t>(chol.n)}, 
                    nn{n * n};
        double *wk_mem1{T::tape->getWKMem(nn)},
               *wk_mem2{T::tape->getWKMem(nn)};
               
        // compute the retun value
        {
            I2 Bi = B;
            for(size_t i = 0; i < nn; ++i, ++Bi)
                wk_mem1[i] = *Bi;
        }
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem1 + i * n); // A^{-1}B
        
        double res_val{0};
        {
            double *diag = wk_mem1;
            for(size_t i = 0; i < n; ++i, diag += n + 1)
                res_val += *diag;
        }
        
        // compute and set the partial derivatives
        for(size_t j = 0; j < n; ++j) // transpose
            for(size_t i = 0; i < n; ++i)
                wk_mem2[i + j * n] = wk_mem1[j + i * n];
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem2 + i * n); // A^{-1}BA^{-1}
        
        T res;
        res.myValue = res_val;
        res.createNode(nn);
        for(size_t i = 0; i < nn; ++i, ++A){
            res.setpDerivatives(i, -wk_mem2[i]);
            res.setpAdjPtrs(i, *A);
        }
        
        return res;
    }
    
    /**
     * same as trInvMatMat_first where the second argument i second iterators is
     * for Ts but the first is not */
    template<class I1, class I2>
    static T trInvMatMat_second(I1 A, I2 B, const CholFactorization &chol){
        static_assert(!is_it_value_type<I1, T>::value,
                      "First iterator is to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        
        const size_t n{static_cast<size_t>(chol.n)}, 
                    nn{n * n};
        double *wk_mem{T::tape->getWKMem(nn)};
               
        // compute the retun value
        {
            I2 Bi = B;
            for(size_t i = 0; i < nn; ++i, ++Bi)
                wk_mem[i] = Bi->value();
        }
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem + i * n); // A^{-1}B
        
        double res_val{0};
        {
            double *diag = wk_mem;
            for(size_t i = 0; i < n; ++i, diag += n + 1)
                res_val += *diag;
        }
        
        // set the partial derivatives
        T res;
        res.myValue = res_val;
        res.createNode(nn);
        double const *A_inv{chol.get_inv()};
        using it_diff_type = typename std::iterator_traits<I2>::difference_type;
        for(size_t j = 0; j < n; ++j, ++A_inv){
            I2 Bi{B + static_cast<it_diff_type>(j * n)};
            I2 Bj{B + static_cast<it_diff_type>(j)};
            for(size_t i = 0; i < j; ++i, ++A_inv, ++Bi, Bj += it_diff_type(n)){
                res.setpDerivatives(i + j * n, *A_inv);
                res.setpDerivatives(j + i * n, *A_inv);
                res.setpAdjPtrs(i + j * n, *Bi);
                res.setpAdjPtrs(j + i * n, *Bj);
            }
            res.setpDerivatives(j * (n + 1), *A_inv);
            res.setpAdjPtrs(j * (n + 1), *Bi);
        }
        
        return res;
    }
    
    /// same as trInvMatMat_first but where both iterators are to Ts
    template<class I1, class I2>
    static T trInvMatMat_identical(I1 A, I2 B, const CholFactorization &chol){
        static_assert(is_it_value_type<I1, T>::value,
                      "First iterator is not to Ts");
        static_assert(is_it_value_type<I2, T>::value,
                      "Second iterator is not to Ts");
        
        const size_t n{static_cast<size_t>(chol.n)}, 
                    nn{n * n};
        double *wk_mem1{T::tape->getWKMem(nn)},
               *wk_mem2{T::tape->getWKMem(nn)};
               
        // compute the retun value
        {
            I2 Bi = B;
            for(size_t i = 0; i < nn; ++i, ++Bi)
                wk_mem1[i] = Bi->value();
        }
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem1 + i * n); // A^{-1}B
        
        double res_val{0};
        {
            double *diag = wk_mem1;
            for(size_t i = 0; i < n; ++i, diag += n + 1)
                res_val += *diag;
        }
        
        // compute and set the partial derivatives
        for(size_t j = 0; j < n; ++j) // transpose
            for(size_t i = 0; i < n; ++i)
                wk_mem2[i + j * n] = wk_mem1[j + i * n];
        for(size_t i = 0; i < n; ++i)
            chol.solve(wk_mem2 + i * n); // A^{-1}BA^{-1}
        
        T res;
        res.myValue = res_val;
        res.createNode(2 * nn);
        double const *A_inv{chol.get_inv()};
        using it_diff_type = typename std::iterator_traits<I2>::difference_type;
        for(size_t j = 0; j < n; ++j, ++A_inv){
            I2 Bi{B + static_cast<it_diff_type>(j * n)};
            I2 Bj{B + static_cast<it_diff_type>(j)};
            for(size_t i = 0; i < j; ++i, ++A_inv, ++Bi, Bj += it_diff_type(n)){
                res.setpDerivatives(i + j * n, *A_inv);
                res.setpDerivatives(j + i * n, *A_inv);
                res.setpAdjPtrs(i + j * n, *Bi);
                res.setpAdjPtrs(j + i * n, *Bj);
            }
            res.setpDerivatives(j * (n + 1), *A_inv);
            res.setpAdjPtrs(j * (n + 1), *Bi);
        }
        
        for(size_t i = 0; i < nn; ++i, ++A){
            res.setpDerivatives(i + nn, -wk_mem2[i]);
            res.setpAdjPtrs(i + nn, *A);
        }
        
        return res;
    }
    
    /**
     * computes the log determinant of a symmetric positive definte matrix
     * using a pre-computed Choleksy decomposition */
    template<class I>
    static T log_deter(I begin, const CholFactorization &chol){
        static_assert(is_it_value_type<I, T>::value,
                      "Iterator is not to Ts");
            
        const size_t n{static_cast<size_t>(chol.n)};
        T res;
        res.myValue = log(chol.determinant());
        res.createNode(n * n);
        
        double const *M_inv{chol.get_inv()};
        using it_diff_type = typename std::iterator_traits<I>::difference_type;
        for(size_t j = 0; j < n; ++j, ++M_inv){
            I Mi{begin + static_cast<it_diff_type>(j * n)};
            I Mj{begin + static_cast<it_diff_type>(j)};
            for(size_t i = 0; i < j; ++i, ++M_inv, ++Mi, Mj += it_diff_type(n)){
                res.setpDerivatives(i + j * n, *M_inv);
                res.setpDerivatives(j + i * n, *M_inv);
                res.setpAdjPtrs(i + j * n, *Mi);
                res.setpAdjPtrs(j + i * n, *Mj);
            }
            res.setpDerivatives(j * (n + 1), *M_inv);
            res.setpAdjPtrs(j * (n + 1), *Mi);
        }
        
        return res;
    }

#endif // if AADLAPACK
};

} // namespace cfaad
