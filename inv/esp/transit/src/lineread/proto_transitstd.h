#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* ../transit/transitstd.c */
extern inline void transitdot P_((int thislevel, int verblevel, ...));
extern int transiterror P_((int flags, const char *str, ...));
extern int vtransiterror P_((int flags, const char *str, va_list ap));
extern int fileexistopen P_((char *in, FILE **fp));
extern int verbfileopen P_((char *in, FILE **fp, char *desc));
extern void transitcheckcalled P_((const long pi, const char *fcn, const int n, ...));
extern void error P_((int exitstatus, int something, const char *fmt, ...));
extern void free_isov P_((prop_isov *isov));
extern void free_isof P_((prop_isof *isof));
extern void free_db P_((prop_db *db));
extern void free_dbnoext P_((prop_dbnoext *db));
extern void free_samp P_((prop_samp *samp));
extern void free_atm P_((prop_atm *atm));
extern void savestr P_((FILE *out, char *str));
extern int reststr P_((FILE *in, char **str));

#undef P_
