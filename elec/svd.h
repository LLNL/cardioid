enum svdMode {svdConverged, svdNotConverged}; 
#ifdef __cplusplus
extern "C" 
{
#endif
int svd(int m, int n, double **a, double *sigma, double **u, double **v);
void svdTest(int m, int n, double **a, double *sigma, double **u, double **v);
void svdLinearLeastSquares(int m, int n, double **a, double *y, double *x) ;
#ifdef __cplusplus
}
#endif
