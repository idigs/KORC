#ifndef NIMROD_SOURCE_H
#define NIMROD_SOURCE_H

#include "fusion_io_source.h"

#ifdef __cplusplus
extern "C" {
#endif
int nimrod_fio_init(const char*);
int nimrod_fio_dealloc();
#ifdef __cplusplus
}
#endif

struct nimrod_search_hint {
  double ixy[2] = {-1.0, -1.0};
  int iblk = 0;
};

class nimrod_source : public fio_source {
public:
  int open(const char*);
  int close();

  int get_available_fields(fio_field_list*) const;
  int get_field_options(fio_option_list*) const;
  int get_field(const field_type, fio_field**, const fio_option_list*);

  int sizeof_search_hint() const
  { return sizeof(nimrod_search_hint); }
  int allocate_search_hint(void** s);
  int deallocate_search_hint(void** s);
};

#endif // NIMROD_SOURCE_H
