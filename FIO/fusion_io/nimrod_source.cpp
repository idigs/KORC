#include <cstring>

#include "fusion_io.h"

int nimrod_source::open(const char* filename)
{
  return nimrod_fio_init(filename);
}

int nimrod_source::close()
{
  return nimrod_fio_dealloc();
}

int nimrod_source::get_available_fields(fio_field_list* fields) const
{
  fields->clear();
  fields->push_back(FIO_DENSITY);
  fields->push_back(FIO_TEMPERATURE);
  fields->push_back(FIO_MAGNETIC_FIELD);
  fields->push_back(FIO_ELECTRIC_FIELD);

  return FIO_SUCCESS;
}

int nimrod_source::get_field_options(fio_option_list* opt) const
{
  opt->clear();
  opt->add_option(FIO_SPECIES, 0);

  return FIO_SUCCESS;
}

int nimrod_source::get_field(const field_type t, fio_field** f,
                             const fio_option_list* opt)
{
  int ierr = -1;
  nimrod_scalar_field* nsf = nullptr;
  switch(t) {
  case FIO_DENSITY:
    nsf = new nimrod_density_field();
    ierr = nsf->load(opt);
    *f = nsf;
    return ierr;
  // case FIO_TEMPERATURE:
  //   nsf = new nimrod_temperature_field();
  //   ierr = nsf->load(opt);
  //   *f = nsf;
  //   return ierr;
  case FIO_MAGNETIC_FIELD:
    *f = new nimrod_magnetic_field();
    return FIO_SUCCESS;
  case FIO_ELECTRIC_FIELD:
    *f = new nimrod_electric_field();
    return FIO_SUCCESS;
  default:
    return FIO_UNSUPPORTED;
  }
}

int nimrod_source::allocate_search_hint(void** s)
{
  *s = new nimrod_search_hint;
  return FIO_SUCCESS;
}

int nimrod_source::deallocate_search_hint(void** s)
{
  delete (nimrod_search_hint*)(*s);
  return FIO_SUCCESS;
}
