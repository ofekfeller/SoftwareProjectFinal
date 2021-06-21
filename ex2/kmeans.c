#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


struct linked_list
{
    double data;
    struct linked_list *next;

};

typedef struct linked_list ELEMENT;
typedef ELEMENT* LINK;

struct points_linked_list
{
    double *vec;
    struct points_linked_list *next;

};

typedef struct points_linked_list MAIN_ELEMENT;
typedef MAIN_ELEMENT* MAIN_LINK;


typedef struct cen_info
{
double *sum;
int cnt;
} CEN_INFO;

typedef CEN_INFO* CEN_LINK;

CEN_LINK init_clusters(int k, int dim){
    CEN_LINK clusters;
    int i;

    clusters = (CEN_INFO*)calloc(k, sizeof(CEN_INFO));
    for ( i = 0; i < k; i++)
    {
        clusters[i].cnt = 0;
        clusters[i].sum = calloc(dim, sizeof(double));
        assert(clusters[i].sum != NULL);
    }

    return clusters;
}
/*
void print_array(double *array, int len){
    int i;

    for ( i = 0; i < len; i++)
    {
        printf("%lf,\n",array[i]);
    }
}
*/
void copy_array(double *source, double *new, int dim){
    int i;
    for (i = 0; i < dim; i++)
    {
        source[i] = new[i];
    }

}

void add_array(double *source, double *new, int dim){
    int i;
    double tmp;

    for (i = 0; i < dim; i++)
    {
        tmp = source[i] + new[i];
        source[i] = tmp;
    }

}

void sub_array(double *source, double *new, int dim){
    int i;
    double tmp;

    for (i = 0; i < dim; i++)
    {
        tmp = source[i] - new[i];
        source[i] = tmp;
    }

}

void update_center(double* origin, CEN_INFO cluster, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        origin[i] = (cluster.sum[i] / cluster.cnt);
    }

}

double norm(double *array1,double *array2, int size){
    double sum;
    double res;
    int i;
    sum=0;

    for (i = 0; i < size; i++)
    {
        res = array1[i]-array2[i];
        sum += (res*res);
    }

    return sum;
}


int classify(double *vec,int num_cent, int size_vec, double **centers){
    int min_ind=0;
    double tmp_norm=0;
    double min_norm;
    int i;

    min_norm = norm(vec, centers[0], size_vec);


    for (i=1; i<num_cent; i++){

        tmp_norm = norm(vec,centers[i],size_vec);

        if (tmp_norm<min_norm){
            min_norm=tmp_norm;
            min_ind=i;
        }
    }
    return min_ind;
}

int update_row(double * to, double * from, int k){
     int i;
    int cnt=0;
    for (i=0; i<k;i++){
        if (from[i]!=to[i]){
            to[i]=from[i];
            cnt++;
        }
    }
    if (cnt>0){
        return 1;
    }
    else{
        return 0;
    }
}

static double** kmeans(double** points, double** centers, int vec_cnt, int counter, int k, int max_iter){

    int i;
    int j;
    int center;
    int *vec_to_cen;
    CEN_LINK clusters;
    int bool;
    int old_center;


    vec_to_cen= calloc(vec_cnt,sizeof(int));
    assert(vec_to_cen != NULL);

    clusters = init_clusters(k, counter);

    for (i = 0; i < vec_cnt; i++)
    {
       center = classify(points[i],  k, counter, centers);

       vec_to_cen[i] = center;

       add_array(clusters[center].sum , points[i], counter);

       clusters[center].cnt++;

    }

    bool = 1;
    for (i=0; i<max_iter-1;i++){

        for(j=0;j<vec_cnt;j++){

            center = classify(points[j] ,k, counter, centers);

            old_center = vec_to_cen[j];

            if (old_center!=center){

                sub_array(clusters[old_center].sum, points[j], counter);

                add_array(clusters[center].sum , points[j], counter);

                clusters[center].cnt++;
                clusters[old_center].cnt--;

                vec_to_cen[j] = center;

                bool = 1;
            }

        }

        if (bool == 0){

            break;

        }

        for (j = 0; j < k; j++)
        {

            update_center(centers[j], clusters[j], counter);

        }

        bool = 0;
    }

    free(vec_to_cen);
    free(clusters);
    free(points);

    return centers;

}

/*
 * API functions
 */


/*
 * Print list of lists of ints without changing it
 */
static PyObject* fit(PyObject *self, PyObject *args)
{

    PyObject *_points, *_centers, *max_iter, *val;
    Py_ssize_t dim, vec_num, i, j, _K, _max_iter;

    if(!PyArg_ParseTuple(args, "OOO", &_points,&_centers, &max_iter)) {
        return NULL;
    }
    if(!PyList_Check(_points) || !PyList_Check(_centers) || !PyList_Check(max_iter)){
        return NULL;
    }

    /* Get the size of it and build the output list*/
    vec_num = PyList_Size(_points);
    dim = PyList_Size(PyList_GetItem(_points, 0));


    double** points = calloc(vec_num, sizeof(double*));
    assert(points!=NULL);
    for(i=0; i<vec_num; i++){
            points[i] = calloc(dim, sizeof(double));
            assert(points[i]!=NULL);
        for(j=0; j<dim; j++){
            val = PyList_GetItem(PyList_GetItem(_points,i), j);
            points[i][j] = PyFloat_AsDouble(val);
        }
    }

    _K = PyList_Size(_centers);
    double** centers = calloc(_K, sizeof(double*));
    assert(centers!=NULL);
    for(i=0; i<_K; i++){
        centers[i] = calloc(dim, sizeof(double));
        assert(centers[i]!=NULL);
        for(j=0; j<dim; j++){
            val = PyList_GetItem(PyList_GetItem(_centers,i), j);
            centers[i][j] = PyFloat_AsDouble(val);
        }
    }

    _max_iter = PyList_Size(max_iter);

    double **clusters = kmeans(points, centers, vec_num,dim, _K, _max_iter);

    PyObject * final_centroids = PyList_New(_K);
    PyObject * centroid;
    for(i=0; i<_K; i++){
        centroid = PyList_New(dim);
        for(j=0; j<dim; j++){

            PyList_SetItem(centroid, j, PyFloat_FromDouble(clusters[i][j]));
        }
        PyList_SetItem(final_centroids, i, centroid);
    }

    free(clusters);
    
    return final_centroids;

}

/*
 * A macro to help us with defining the methods
 * Compare with: {"f1", (PyCFunction)f1, METH_NOARGS, PyDoc_STR("No input parameters")}
*/
#define FUNC(_flag, _name, _docstring) { #_name, (PyCFunction)_name, _flag, PyDoc_STR(_docstring) }

static PyMethodDef _methods[] = {
    FUNC(METH_VARARGS, fit, "Find k means clusters given points, early centers and some shitty integers"),
    {NULL, NULL, 0, NULL}   /* sentinel */
};

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    return PyModule_Create(&_moduledef);
}


