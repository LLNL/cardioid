#ifndef CS_GETTIME__
#define CS_GETTIME__

#ifdef __cplusplus
extern "C" {
#endif

  void cs_inittime(void);
  double cs_gettime(void);
  double cs_getlocaloffset(void);

#ifdef __cplusplus
extern "C" {
#endif

#endif
