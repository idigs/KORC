#ifndef NIMROD_FIELD_H
#define NIMROD_FIELD_H

#include "fusion_io_field.h"
#include "nimrod_source.h"
#include "options.h"

#ifdef __cplusplus
extern "C" {
#endif
int nimrod_fio_get_time(const double*);
int nimrod_fio_ndens_eval(const field_attribute, const double*, double*);
int nimrod_fio_ndens_eval_hint(const field_attribute, const double*, double*,
                                nimrod_search_hint*);
// int nimrod_fio_temp_eval(const field_attribute, const double*, double*);
// int nimrod_fio_temp_eval_hint(const field_attribute, const double*, double*,
//                                nimrod_search_hint*);
int nimrod_fio_b_eval(const double*, double*);
int nimrod_fio_b_eval_hint(const double*, double*, nimrod_search_hint*);
// int nimrod_fio_b_eval_deriv(const double*, double*);
// int nimrod_fio_b_eval_deriv_hint(const double*, double*, nimrod_search_hint*);
int nimrod_fio_e_eval(const double*, double*);
int nimrod_fio_e_eval_hint(const double*, double*, nimrod_search_hint*);
// int nimrod_fio_e_eval_deriv(const double*, double*);
// int nimrod_fio_e_eval_deriv_hint(const double*, double*, nimrod_search_hint*);
#ifdef __cplusplus
}
#endif

class nimrod_fio_field : public fio_field {
 public:
  int get_real_parameter(const field_parameter, double*);
};

class nimrod_scalar_field : public nimrod_fio_field {
 protected:
  field_attribute species = 0;
 public:
  int dimension() const { return 1; }
  int load(const fio_option_list*);
};

class nimrod_vector_field : public nimrod_fio_field {
 public:
  int dimension() const { return 3; }
};

class nimrod_density_field : public nimrod_scalar_field {
 public:
  virtual int eval(const double*, double*, void* = nullptr);
  virtual fio_field* clone() const { return new nimrod_density_field(*this); }
};

// class nimrod_temperature_field : public nimrod_scalar_field {
//  public:
//   virtual int eval(const double*, double*, void* = nullptr);
//   virtual fio_field* clone() const { return new nimrod_temperature_field(*this); }
// };

class nimrod_magnetic_field : public nimrod_vector_field {
 public:
  virtual int eval(const double*, double*, void* = nullptr);
  // virtual int eval_deriv(const double*, double*, void* = nullptr);
  virtual fio_field* clone() const { return new nimrod_magnetic_field(*this); }
};

class nimrod_electric_field : public nimrod_vector_field {
 public:
  virtual int eval(const double*, double*, void* = nullptr);
  // virtual int eval_deriv(const double*, double*, void* = nullptr);
  virtual fio_field* clone() const { return new nimrod_electric_field(*this); }
};

#endif // NIMROD_FIELD_H
