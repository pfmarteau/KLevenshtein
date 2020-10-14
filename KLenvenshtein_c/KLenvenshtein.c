/*  Wrapping cos function from math.h with the Python-C-API. */

#include <Python.h>
#include <math.h>

static double **localKernel(const char *A, int la, const char *B, int lb, double sigma, double epsilon){
    double **LK = (double **)calloc(la, sizeof(double *));
    double e0 = 1.0/3.0;
    double e1 = (exp(-1.0/sigma)+epsilon)/(3*(1+epsilon));
    for (int i = 0; i < la; i++) {
    	LK[i] = (double *)calloc(lb, sizeof(double));
        for (int j=0; j<lb; j++) {
                if (A[i] == B[j])
			LK[i][j] = e0;
                else
			LK[i][j] = e1;
        }
    }
    return LK;
}

static double **localKernel_int(const int *A, int la, const int *B, int lb, double sigma, double epsilon){
    double **LK = (double **)calloc(la, sizeof(double *));
    double e0 = 1.0/3.0;
    double e1 = (exp(-1.0/sigma)+epsilon)/(3.0*(1+epsilon));
    for (int i = 0; i < la; i++) {
    	LK[i] = (double *)calloc(lb, sizeof(double));
        for (int j=0; j<lb; j++) {
                if (A[i] == B[j])
			LK[i][j] = e0;
                else
			LK[i][j] = e1;
        }
    }
    return LK;
}

double kLevenshtein_lk(int lA, int lB, double **local_kernel){
    int la = lA+1;
    int lb = lB+1;
    int lmx = la;
    int l = lb;
    if (lb > la){
	lmx = lb;
        l = la;
    }
    double **DP = (double **)calloc(la, sizeof(double *));
    double **DP1 = (double **)calloc(la, sizeof(double *));
    for (int i = 0; i < la; i++) {
    	DP[i] = (double *)calloc(lb, sizeof(double));
	DP1[i] = (double *)calloc(lb, sizeof(double));
    }

    double *DP2 = (double *)calloc(lmx, sizeof(double));

    DP2[0]=1.0;
    for (int i=1; i<l; i++)
        DP2[i] = local_kernel[i-1][i-1];

    DP[0][0] = 1;
    DP1[0][0] = 1;

    for (int i=1; i<la; i++){
        DP[i][1] = DP[i-1][1]*local_kernel[i-1][2];
        DP1[i][1] = DP1[i-1][1]*DP2[i];
        }
    for (int j=1; j<lb; j++){
        DP[1][j] = DP[1][j-1]*local_kernel[2][j-1];
        DP1[1][j] = DP1[1][j-1]*DP2[j];
        }

    for (int i = 1; i<la; i++){
        for (int j=1; j<lb; j++){ 
            DP[i][j] = (DP[i-1][j] + DP[i][j-1] + DP[i-1][j-1]) * local_kernel[i-1][j-1];
            if (i == j)
                DP1[i][j] = DP1[i-1][j-1] * local_kernel[i-1][j-1] + DP1[i-1][j] * DP2[i] + DP1[i][j-1]  *DP2[j];
            else
                DP1[i][j] = DP1[i-1][j] * DP2[i] + DP1[i][j-1] * DP2[j];
    	}
    }
    double ans=DP[la-1][lb-1]+DP1[la-1][lb-1];
    for (int i=0; i<la; i++){
	free(DP[i]);
        free(DP1[i]);
    }
    free(DP);
    free(DP1);
    free(DP2);
    return ans;
}

int *getIntArray(PyObject* seq){
    /* prepare data as an array of doubles */
    int seqlen = PySequence_Fast_GET_SIZE(seq);
    int *iarr = malloc(seqlen*sizeof(int));
    if(!iarr) {
        Py_DECREF(seq);
        return NULL;
    }
    for(int i=0; i < seqlen; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
        if(!item) {
            Py_DECREF(seq);
            free(iarr);
            return NULL;
        }
        fitem = PyNumber_Float(item);
        if(!fitem) {
            Py_DECREF(seq);
            free(iarr);
            PyErr_SetString(PyExc_TypeError, "all items must be numbers");
            return 0;
        }
        iarr[i] = PyLong_AsLong(fitem);
        Py_DECREF(fitem);
    }    
    return iarr;
}

static PyObject* similarity_str(PyObject* self, PyObject* args)
{
    double answer;
    char *A;
    char *B;
    double sigma, epsilon;
    if (!PyArg_ParseTuple(args, "ssdd", &A, &B, &sigma, &epsilon)) {
        printf("wrong args: expect string1, string2, (float) sigma, (float) epsilon \n");
    return Py_BuildValue("f", -1.0);
    }

    int argc = PyTuple_GET_SIZE(args);
    if (argc!=4){
        printf("wrong number of args, string1, string2, (float) sigma, (float) epsilon \n");
        return Py_BuildValue("f", -1.0);
        }

    int la = strlen(A);
    int lb = strlen(B);
    double **LK = localKernel(A, la, B, lb, sigma, epsilon);
    
    answer = kLevenshtein_lk(la, lb, LK);
    for (int i=0; i<la; i++)
       free(LK[i]);
    free(LK);
    
    /*  construct the output from cos, from c double to python float */
    return Py_BuildValue("f", answer);
}

static PyObject* similarity_int(PyObject* self, PyObject* args)
{
    double answer;
    int *A;
    int *B;
    double sigma, epsilon;

    int argc = PyTuple_GET_SIZE(args);

    PyObject *seq1, *seq2;
    if(!PyArg_ParseTuple(args, "OOdd", &seq1, &seq2, &sigma, &epsilon)){
 	printf("wrong number of args, string1, string2, (float) sigma, (float) epsilon \n");
    return Py_BuildValue("f", -1.0);
    }

    if (argc!=4)
        printf("wrong number of args, string1, string2, (float) sigma, (float) epsilon \n");

    A = getIntArray(seq1);
    B = getIntArray(seq2);

    if (argc!=4){
        printf("wrong number of args, string1, string2, (float) sigma, (float) epsilon \n");
        return Py_BuildValue("f", -1.0);
        }

   int la = PySequence_Fast_GET_SIZE(seq1);
   int lb = PySequence_Fast_GET_SIZE(seq2);
   double **LK = localKernel_int(A, la, B, lb, sigma, epsilon);
   free(A);
   free(B);
   answer = kLevenshtein_lk(la, lb, LK);
   for (int i=0; i<la; i++)
       free(LK[i]);
   free(LK);
   /*  construct the output, from c double to python float */
   //printf("---%f \n", answer);
   return Py_BuildValue("d", answer);
}

/*  define functions in module */
static PyMethodDef similarityMethods[] =
{
     {"similarity_str", similarity_str, METH_VARARGS, "evaluate kLevenshtein with string args"},
     {"similarity_int", similarity_int, METH_VARARGS, "evaluate kLevenshtein with array of int args"},
     {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3
/* module initialization */
/* Python version 3*/
static struct PyModuleDef cModPyDem =
{
    PyModuleDef_HEAD_INIT,
    "KLevenshtein", "Some documentation",
    -1,
    similarityMethods
};

PyMODINIT_FUNC
PyInit_KLevenshtein(void)
{
    return PyModule_Create(&cModPyDem);
}

#else

/* module initialization */
/* Python version 2 */
PyMODINIT_FUNC
initKLevenshtein(void)
{
    (void) Py_InitModule("KLevenshtein", similarityMethods);
}


#endif
