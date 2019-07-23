#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void c_dged(void *, void *, void *, void *, void *, void *, void *);
extern void c_dgh(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dghst(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dghyp(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_djsu(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dsged(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dsnig(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dsnorm(void *, void *, void *, void *, void *, void *, void *);
extern void c_dsstd(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_dstd(void *, void *, void *, void *, void *, void *, void *);
extern void c_pged(void *, void *, void *, void *, void *, void *);
extern void c_pjsu(void *, void *, void *, void *, void *, void *, void *);
extern void c_psged(void *, void *, void *, void *, void *, void *, void *);
extern void c_psnorm(void *, void *, void *, void *, void *, void *);
extern void c_psstd(void *, void *, void *, void *, void *, void *, void *);
extern void c_pstd(void *, void *, void *, void *, void *, void *);
extern void c_qged(void *, void *, void *, void *, void *, void *);
extern void c_qjsu(void *, void *, void *, void *, void *, void *, void *);
extern void c_qsged(void *, void *, void *, void *, void *, void *, void *);
extern void c_qsnorm(void *, void *, void *, void *, void *, void *);
extern void c_qsstd(void *, void *, void *, void *, void *, void *, void *);
extern void c_qstd(void *, void *, void *, void *, void *, void *);
extern void c_rged(void *, void *, void *, void *, void *);
extern void c_rghst(void *, void *, void *, void *, void *, void *);
extern void c_rghyp(void *, void *, void *, void *, void *, void *, void *);
extern void c_rjsu(void *, void *, void *, void *, void *, void *);
extern void c_rsged(void *, void *, void *, void *, void *, void *);
extern void c_rsnig(void *, void *, void *, void *, void *, void *);
extern void c_rsnorm(void *, void *, void *, void *, void *);
extern void c_rsstd(void *, void *, void *, void *, void *, void *);
extern void c_rstd(void *, void *, void *, void *, void *);
extern void garchfilterC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void garchsimC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void mixturefilterC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void qNIG(void *, void *, void *, void *, void *, void *, void *);
extern void rfilter(void *, void *, void *, void *);
extern void star3dfun(void *, void *, void *, void *);
extern void starxfilters1s(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters2s(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters2sx(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters3s(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters3sx(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters4s(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxfilters4sx(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxsim1(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxsim2(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxsim3(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void starxsim4(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
  {"c_dged",          (DL_FUNC) &c_dged,           7},
  {"c_dgh",           (DL_FUNC) &c_dgh,            9},
  {"c_dghst",         (DL_FUNC) &c_dghst,          8},
  {"c_dghyp",         (DL_FUNC) &c_dghyp,          9},
  {"c_djsu",          (DL_FUNC) &c_djsu,           8},
  {"c_dsged",         (DL_FUNC) &c_dsged,          8},
  {"c_dsnig",         (DL_FUNC) &c_dsnig,          8},
  {"c_dsnorm",        (DL_FUNC) &c_dsnorm,         7},
  {"c_dsstd",         (DL_FUNC) &c_dsstd,          8},
  {"c_dstd",          (DL_FUNC) &c_dstd,           7},
  {"c_pged",          (DL_FUNC) &c_pged,           6},
  {"c_pjsu",          (DL_FUNC) &c_pjsu,           7},
  {"c_psged",         (DL_FUNC) &c_psged,          7},
  {"c_psnorm",        (DL_FUNC) &c_psnorm,         6},
  {"c_psstd",         (DL_FUNC) &c_psstd,          7},
  {"c_pstd",          (DL_FUNC) &c_pstd,           6},
  {"c_qged",          (DL_FUNC) &c_qged,           6},
  {"c_qjsu",          (DL_FUNC) &c_qjsu,           7},
  {"c_qsged",         (DL_FUNC) &c_qsged,          7},
  {"c_qsnorm",        (DL_FUNC) &c_qsnorm,         6},
  {"c_qsstd",         (DL_FUNC) &c_qsstd,          7},
  {"c_qstd",          (DL_FUNC) &c_qstd,           6},
  {"c_rged",          (DL_FUNC) &c_rged,           5},
  {"c_rghst",         (DL_FUNC) &c_rghst,          6},
  {"c_rghyp",         (DL_FUNC) &c_rghyp,          7},
  {"c_rjsu",          (DL_FUNC) &c_rjsu,           6},
  {"c_rsged",         (DL_FUNC) &c_rsged,          6},
  {"c_rsnig",         (DL_FUNC) &c_rsnig,          6},
  {"c_rsnorm",        (DL_FUNC) &c_rsnorm,         5},
  {"c_rsstd",         (DL_FUNC) &c_rsstd,          6},
  {"c_rstd",          (DL_FUNC) &c_rstd,           5},
  {"garchfilterC",    (DL_FUNC) &garchfilterC,    16},
  {"garchsimC",       (DL_FUNC) &garchsimC,       12},
  {"mixturefilterC",  (DL_FUNC) &mixturefilterC,  10},
  {"qNIG",            (DL_FUNC) &qNIG,             7},
  {"rfilter",         (DL_FUNC) &rfilter,          4},
  {"star3dfun",       (DL_FUNC) &star3dfun,        4},
  {"starxfilters1s",  (DL_FUNC) &starxfilters1s,  14},
  {"starxfilters2s",  (DL_FUNC) &starxfilters2s,  14},
  {"starxfilters2sx", (DL_FUNC) &starxfilters2sx, 12},
  {"starxfilters3s",  (DL_FUNC) &starxfilters3s,  14},
  {"starxfilters3sx", (DL_FUNC) &starxfilters3sx, 12},
  {"starxfilters4s",  (DL_FUNC) &starxfilters4s,  14},
  {"starxfilters4sx", (DL_FUNC) &starxfilters4sx, 12},
  {"starxsim1",       (DL_FUNC) &starxsim1,       10},
  {"starxsim2",       (DL_FUNC) &starxsim2,       10},
  {"starxsim3",       (DL_FUNC) &starxsim3,       10},
  {"starxsim4",       (DL_FUNC) &starxsim4,       10},
  {NULL, NULL, 0}
};

void R_init_twinkle(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
