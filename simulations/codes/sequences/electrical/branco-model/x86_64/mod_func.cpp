#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _ar_reg(void);
extern void _cad_reg(void);
extern void _caL3d_reg(void);
extern void _ca_reg(void);
extern void _CaT_reg(void);
extern void _HH2_reg(void);
extern void _h_reg(void);
extern void _inwardrect_reg(void);
extern void _kca_reg(void);
extern void _kir_reg(void);
extern void _km_reg(void);
extern void _kv_reg(void);
extern void _na_reg(void);
extern void _NMDA_Mg_T_reg(void);
extern void _release_BMK_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"mod.files/ar.mod\"");
    fprintf(stderr, " \"mod.files/cad.mod\"");
    fprintf(stderr, " \"mod.files/caL3d.mod\"");
    fprintf(stderr, " \"mod.files/ca.mod\"");
    fprintf(stderr, " \"mod.files/CaT.mod\"");
    fprintf(stderr, " \"mod.files/HH2.mod\"");
    fprintf(stderr, " \"mod.files/h.mod\"");
    fprintf(stderr, " \"mod.files/inwardrect.mod\"");
    fprintf(stderr, " \"mod.files/kca.mod\"");
    fprintf(stderr, " \"mod.files/kir.mod\"");
    fprintf(stderr, " \"mod.files/km.mod\"");
    fprintf(stderr, " \"mod.files/kv.mod\"");
    fprintf(stderr, " \"mod.files/na.mod\"");
    fprintf(stderr, " \"mod.files/NMDA_Mg_T.mod\"");
    fprintf(stderr, " \"mod.files/release_BMK.mod\"");
    fprintf(stderr, "\n");
  }
  _ar_reg();
  _cad_reg();
  _caL3d_reg();
  _ca_reg();
  _CaT_reg();
  _HH2_reg();
  _h_reg();
  _inwardrect_reg();
  _kca_reg();
  _kir_reg();
  _km_reg();
  _kv_reg();
  _na_reg();
  _NMDA_Mg_T_reg();
  _release_BMK_reg();
}

#if defined(__cplusplus)
}
#endif
