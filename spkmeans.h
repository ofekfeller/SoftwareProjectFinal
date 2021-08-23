int classify(double *vec,int num_cent, int size_vec, double** centers);
void copy_array(double *from, double* to, int dim);

typedef struct cen_info
{
    double *sum; 
    int cnt;
} CEN_INFO;
typedef CEN_INFO* CEN_LINK;

CEN_LINK init_clusters(int k, int dim);
void update_center(double* origin, CEN_INFO cluster, int dim);
double* sub_vectors(double *vec1,double *vec2, int n);
void add_vectors(double *vec1, double *vec2, int n);
void print_vec(double* vec, int len);
void print_mat(double** mat, int n, int m);
double norm(double *vec, int size);
double p_exp(double *vec1,double *vec2, int size);
void normalize(double* vec, int size);
void normalize_mat(double** mat, int rows, int cols);
double** matrix_sub(double** mat1, double** mat2, int n);
double** matrix_mul(double** mat1, double** mat2, int x, int y, int z);
double** sq_matrix_mul(double** mat1, double** mat2, int n);
void mul_lines(double** mat, double* vec, int dim);
void mul_columns(double** mat, double* vec, int dim);
double** weighted_matrix(double** points, int dim, int vec_num);
double sum_line(double **mat, int n, int index);
double* vector_sum(double** mat, int n);
void div_square_vec(double* vec, int dim);
double** get_diag_mat(double* vec, int n);
double** transpose(double **mat, int n, int m);
double** wrap(double **mat1,double **mat2, int n);
double* get_diag(double** mat, int n);
double get_abs(double num);
int* get_max_indexes(double** mat, int n);
double compute_theta(double x, double y, double z);
double sign(double x);
double compute_t(double x);
double compute_c(double x);
double square_sum(double** mat, int n);
double compute_off(double** mat, int n);
double** many_mul(double*** mat, int n, int m);
double** get_eye_mat(int dim);
double** create_p_matrix(double s, double c, int x, int y, int n);
double* get_deltas(double* lst, int n);
int determine_k(double* lst, int n);
void swap(int *a, int *b);
int partition(int *array, int low, int high);
void quickSort(int *array, int low, int high);
int* read_file_dimensions(char* filename);
double** init_2d_array(int n, int m);
double** get_points_from_file(char* filename, int vec_len, int vec_num);
double* get_diag_vec(double** weighted, int dim);
double** get_normalized_matrix(double** weighted, double* diag, int dim);
double** deep_copy(double** mat, int dim);
void compute_normalized(double** mat, int dim, double c, double s, int i, int j);

typedef struct eigen_ret{
    double * eigen_values;
    double** eigen_vectors;
    int k;
} EIGEN;
typedef EIGEN* EIGEN_LINK;

EIGEN_LINK get_eigens_and_k(double** normalized, int dim, int k);
void kmeans_goal(double** points, char* goal, int vec_num, int dim);
double** get_spk_points(double** points, int dim, int vec_num, int k);
