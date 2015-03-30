#include <R_ext/Applic.h>
#include <Rcpp.h>
#include <cmath>
#include <stdexcept>
#include <iterator>
#include <functional>
#include <numeric>

#include "optimize.h"

namespace {

  /* represent the variables used in the iteration */

  struct item_data {
    double flag, count, lambda, log_lambda;
    item_data(const double& flag, const double& count,
	      const double& lambda, const double& log_lambda) :
      flag(flag), count(count),
      lambda(lambda), log_lambda(log_lambda) {}
  };

  /* A *much* better version of this class is in BOOST */

  template <class T>
  class four_iterator : public std::iterator<std::input_iterator_tag,
					     item_data> {
    typedef T base_iterator;
    base_iterator flag_it, count_it, lambda_it, log_lambda_it;
  public:
    four_iterator() {}
    four_iterator(const base_iterator& flag_it,
		  const base_iterator& count_it,
		  const base_iterator& lambda_it,
		  const base_iterator& log_lambda_it) :
      flag_it(flag_it), count_it(count_it),
      lambda_it(lambda_it), log_lambda_it(log_lambda_it) {}

    four_iterator(const four_iterator<T>& other) :
      flag_it(other.flag_it), count_it(other.count_it),
      lambda_it(other.lambda_it), log_lambda_it(other.log_lambda_it) {}

    four_iterator& operator++() {
      ++flag_it; ++count_it; ++lambda_it; ++log_lambda_it;
      return *this;
    }

    four_iterator operator++(int) {
      four_iterator tmp(*this);
      operator++();
      return tmp;
    }

    bool operator==(const four_iterator<T>& other) {
      return flag_it == other.flag_it;
    }

    bool operator!=(const four_iterator<T>& other) {
      return flag_it != other.flag_it;
    }

    item_data operator*() {
      return item_data(*flag_it, *count_it, *lambda_it, *log_lambda_it);
    }

  };

  /*
    This is an input iterator that transforms results as it goes.
    Boost (again) does this better.
  */

  template <class T, class it>
  class transform_iterator : public std::iterator<std::input_iterator_tag,
						  typename T::result_type> {
  private:
    typedef it base_iterator_t;
    typedef T transformer_t;
    transformer_t transform;
    base_iterator_t base;
    typedef typename transformer_t::result_type my_res_t;
  public:
    transform_iterator() {}

    transform_iterator(const transformer_t& transform,
		       const base_iterator_t& base) :
      transform(transform), base(base) {}

    transform_iterator& operator++() {
      ++base;
      return *this;
    }

    transform_iterator operator++(int) {
      transform_iterator tmp(*this);
      operator++();
      return tmp;
    }

    bool operator==(const transform_iterator& other) {
      return base == other.base;
    }

    bool operator!=(const transform_iterator& other) {
      return base != other.base;
    }

    my_res_t operator*() {
      return transform(*base);
    }

  };

  /*
    Herewith a functor that computes contributions to the cost function
  */

  class transform_to_cost : public std::unary_function<item_data, double> {
  private:
    double delta, delta_log_delta, lgamma_delta;
    double duration, log_duration;
  public:
    transform_to_cost(const double& delta, const double& duration) :
      delta(delta),
      delta_log_delta(delta * std::log(delta)),
      lgamma_delta(lgamma(delta)),
      duration(duration), log_duration(std::log(duration)) {}

    double operator()(const item_data& datum) const {
      const double shift_count = datum.count + delta;
      const double lcomb = lgamma(shift_count) -
	(lgamma_delta + lgamma(1 + datum.count));
      return (datum.count * (log_duration + datum.log_lambda)
	      + delta_log_delta + lcomb
	      - shift_count * std::log(delta + datum.lambda * duration));
    }

  };

  /*
     here we handle contributions to the second derivative matrix.
     For various technical reasons there are just two independent
     terms, so we only think about those.
  */

  struct hessian_contrib {
    double d2_dintercept2;
    double d2_dgradient2;

    hessian_contrib() : d2_dintercept2(), d2_dgradient2() {}

    hessian_contrib(const double& d2_dintercept2,
		    const double& d2_dgradient2) :
      d2_dintercept2(d2_dintercept2),
      d2_dgradient2(d2_dgradient2) {}

    hessian_contrib& operator+=(const hessian_contrib& other) {
      d2_dintercept2 += other.d2_dintercept2;
      d2_dgradient2 += other.d2_dgradient2;
      return *this;
    }

    hessian_contrib operator+(const hessian_contrib& other) {
      hessian_contrib res(*this);
      res += other;
      return res;
    }
  };

  /* and a functor that computes contributions to the hessian */

  class transform_to_hessian :
    public std::unary_function<item_data, hessian_contrib> {
  private:
    double delta, delta_term_1, delta_term_2;
    double duration, log_duration;
  public:
    transform_to_hessian(const double& delta, const double& duration) :
      delta(delta),
      delta_term_1(std::log(delta) + 1 - R::digamma(delta)),
      delta_term_2(1/delta - R::trigamma(delta)),
      duration(duration), log_duration(std::log(duration)) {}

    hessian_contrib operator()(const item_data& datum) const {
      const double shift_count = datum.count + delta;
      const double shift_delta = delta + datum.lambda * duration;
      const double d_dlambda = datum.count / datum.lambda -
	shift_count * duration / shift_delta;
      const double d_dloglambda = datum.lambda * d_dlambda;
      const double d2_dlambda2 = - datum.count /  (datum.lambda * datum.lambda)
	+ shift_count * duration * duration / (shift_delta * shift_delta);
      const double d2_dloglambda2 = datum.lambda * datum.lambda * d2_dlambda2
	+ d_dloglambda;
      return hessian_contrib(d2_dloglambda2,
			     datum.flag * d2_dloglambda2);
    }

  };

  class CostData : public std::unary_function<double, double> {
  private:
    typedef four_iterator<Rcpp::NumericVector::const_iterator> iterator_type;
  public:
    CostData(const Rcpp::NumericVector count,
	     const Rcpp::NumericVector flag,
	     const Rcpp::NumericVector log_lambda,
	     const double duration) :
      count(count), flag(flag), log_lambda(log_lambda),
      duration(duration),
      lambda(Rcpp::exp(log_lambda)) {
      if (count.size() != flag.size()) {
	throw std::invalid_argument("count and flag mismatch");
      }
      if (log_lambda.size() != flag.size()) {
	throw std::invalid_argument("log_lambda and flag mismatch");
      }
      if (duration <= 0) {
	throw std::invalid_argument("negative duration");
      }
    }

    double operator()(const double epsilon) const {
      const double delta = 1/epsilon;
      if (delta <= 0) {
	return R_PosInf;
      }

      iterator_type begin = iterator_type(flag.begin(),
					  count.begin(),
					  lambda.begin(),
					  log_lambda.begin());
      iterator_type end = iterator_type(flag.end(),
					count.end(),
					lambda.end(),
					log_lambda.end());

      transform_to_cost transformer(delta, duration);
      typedef transform_iterator<transform_to_cost, iterator_type> xform_t;
      return std::accumulate(xform_t(transformer, begin),
			     xform_t(transformer, end),
			     double(),
			     std::minus<double>());
    }


    double best(const double startx) const {
      double start = startx, start_val = operator()(start);
      double low = 0.5 * start, low_val = operator()(low);
      while (low_val < start_val) {
	start = low; start_val = low_val;
	low *= 0.5; low_val = operator()(low);
      }
      double high = 2 * start, high_val = operator()(high);
      while (high_val < start_val) {
	start = high; start_val = high_val;
	high *= 2; high_val = operator()(high);
      }
      return minimize(*this, low, high);
    }

    Rcpp::NumericMatrix covariance(const double epsilon) const {
      const double delta = 1 / epsilon;
      iterator_type begin = iterator_type(flag.begin(),
					  count.begin(),
					  lambda.begin(),
					  log_lambda.begin());
      iterator_type end = iterator_type(flag.end(),
					count.end(),
					lambda.end(),
					log_lambda.end());
      transform_to_hessian transformer(delta, duration);
      typedef transform_iterator<transform_to_hessian, iterator_type> xform_t;
      hessian_contrib summed = std::accumulate(xform_t(transformer, begin),
					       xform_t(transformer, end),
					       hessian_contrib());
      Rcpp::NumericMatrix res(2, 2);
      const double t1 = -summed.d2_dintercept2;
      const double t2 = -summed.d2_dgradient2;
      const double d2 = delta * delta;
      const double det = 1 / (t2 * (t1 - t2));
      res(0, 0) = t2 * det;
      res(1, 1) = t1 * det;
      res(0, 1) = res(1, 0) = -t2 * det;
      return res;
    }
  private:
    Rcpp::NumericVector count, flag, log_lambda, lambda;
    double duration;
  };


  SEXP fitNegBinData(SEXP coeff, SEXP count, SEXP flag,
		     SEXP log_lambda, SEXP duration) {
    BEGIN_RCPP
      Rcpp::NumericVector res(Rcpp::clone(coeff));
    CostData cost(count, flag, log_lambda, Rcpp::as<double>(duration));
    if (res.size() != 1) {
      throw std::invalid_argument("Bad size for coefficients.");
    }

    res[0] = cost.best(res[0]);

    Rcpp::List result;
    result["par"] = res;
    result["cov"] = cost.covariance(res[0]);
    return result;
    END_RCPP
      }
}


static R_CallMethodDef callMethods[] = {
  {".call.fitNegBinData", (DL_FUNC) &fitNegBinData, 5},
  NULL
};

RcppExport void R_init_assurance(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
