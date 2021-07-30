
#include <Math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int epsilon = 0.001;

double* sub_vectors(double *vec1,double *vec2, int n){
    int i;
    double* vec=calloc(n,sizeof(double));
    assert(vec!=NULL);
    for (i=0; i<n;i++){
        vec[i]=vec1[i]-vec2[i];
    }
    return vec;
}

void print_vec(double* vec, int len){
    char sign;
    for(int i=0;i<len;i++){
            if (i==len-1) sign='\n';
            else sign = ',';
            printf("%lf%c", vec[i], sign);
    }
}

void print_mat(double** mat, int n, int m){
    char sign;

    for (int i=0; i<n; i++){
    
        for (int j=0; j < m; j++)
        {
            if (j==m-1) sign='\n';
            else sign = ',';
            
            printf("%lf%c", mat[i][j], sign);
        }
    }

}

/*
compute norm of a vector
*/
double norm(double *vec, int size){
    double sum;
    double res;
    int i;
    sum=0;

    for (i = 0; i < size; i++)
    {
        sum += (vec[i]*vec[i]);
    }
    
    return sqrt(sum);
}

/*
given the norm of two vectors, will computer the weight in the matrix
*/
double p_exp(double *vec1,double *vec2, int size){
    double my_norm;
    double *vec=sub_vectors(vec1,vec2,size);
    my_norm=norm(vec, size);
    my_norm=-my_norm/2;
    return exp(my_norm);
}
/*
divide each cordinate of the vector by the norm of the vector
*/
void normalize(double* vec, int size){
    int i;
    double vec_norm;
    vec_norm=norm(vec,size);
    for (i=0;i<size;i++){
        vec[i]=vec[i]/vec_norm;
    }
}

// normalize matrix by its rows
void normalize_mat(double** mat, int rows, int cols){
    for(int i=0;i<rows;i++){
        normalize(mat[i], cols);
    }
}

/*
return new matrix that is a substraction of two matrixes
*/
double** matrix_sub(double** mat1, double** mat2, int n){
    double *a;
    double **mat;
    int i;
    int j;
    a=calloc(n*n,sizeof(double));
    assert(a!=NULL);
    mat=calloc(n,sizeof(double*));
    assert(mat!=NULL);
    for (i=0;i<n;i++){
        mat[i]= a + i*n;
    }
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            mat[i][j]= mat1[i][j] - mat2[i][j];
        }
    }
    return mat;
}

/*
return new matrix that is a multiplication of two matrixes
*/
double** matrix_mul(double** mat1, double** mat2, int x, int y, int z){
    double *a;
    double **mat;
    int i,j,k;
    a=calloc(x*z , sizeof(double));
    assert(a!=NULL);
    mat=calloc(x,sizeof(double*));
    assert(mat!=NULL);
    for (i=0;i<x;i++){
        mat[i]= a + i*z;
    }

    for (i = 0; i < x; i++) {
        for (j = 0; j < z; j++) {
            mat[i][j] = 0;
            for (k = 0; k < y; k++)
                mat[i][j] += mat1[i][k] * mat2[k][j];
        }
    }
    return mat;
}
/*
return new matrix that is a multiplication of two square matrixes
*/
double** sq_matrix_mul(double** mat1, double** mat2, int n){
    return matrix_mul(mat1, mat2, n, n, n);
}

/*
given a vector and matrix- multiply each matrix row with the corresponding vector value
*/
void mul_lines(double** mat, double* vec, int dim){
    for(int i=0; i<dim; i++){
          for(int j=0; j<dim; j++){
              mat[i][j] *= vec[i];
    }
}
}

/*
given a vector and matrix- multiply each matrix column with the corresponding vector value
*/
void mul_columns(double** mat, double* vec, int dim){
    for(int i=0; i<dim; i++){
          for(int j=0; j<dim; j++){
              mat[i][j] *= vec[j];
    }
}
}
/*
given set of points compute the weight between each and store them in wheightened adjacency matrix
*/
double** weighted_matrix(double** points, int dim, int vec_num){
    double *a;
    double **mat;
    int i;
    int j;
    a=calloc(vec_num*vec_num,sizeof(double));
    assert(a!=NULL);
    mat=calloc(vec_num,sizeof(double*));
    assert(mat!=NULL);
    for (i=0;i<vec_num;i++){
        mat[i]= a + i*vec_num;
    }
    for (i=0; i<vec_num;i++){
        for (j=i; j<vec_num; j++){
            mat[i][j]=p_exp(points[i],points[j],dim);
            mat[j][i]=mat[i][j];
        }
    }
    return mat;
}
/*
given a matrix and an index, return the sum of the elements in the row of the index in the matrix
*/
double sum_line(double **mat, int n, int index){
    int sum=0;
    int i;
    for (i=0; i<n; i++){
        sum+=mat[index][i];
    }
    return sum;
}
/*
given a matrix, return a vector where each cordinate i, is the sum of row i in the matrix
*/
double* vector_sum(double** mat, int n){
    double* vec;
    int i;
    vec=calloc(n,sizeof(double));
    assert(vec!=NULL);
    for (i=0;i<n;i++){
        vec[i]=sum_line(mat,n,i);
    }
    return vec;
}


/*
given vector, calculate for each val - 1/sqrt(val)
*/
void div_square_vec(double* vec, int dim){
    for(int i=0; i<dim; i++){
        vec[i] = 1/sqrt(vec[i]);
    }
}


/*
given a vector of sums, return the diagonal degree matrix
*/
double** get_diag_mat(double* vec, int n){  //yoni please make sure you dont pass the same indexes twice
    double *a;
    double **mat;
    int i;
    int j;
    a=calloc(n*n,sizeof(double));
    assert(a!=NULL);
    mat=calloc(n,sizeof(double*));
    assert(mat!=NULL);
    for (i=0;i<n;i++){
        mat[i]= a + i*n;
    }
    for (i=0; i<n; i++){
        for (j=0; j<n;j++){
            if (i==j){
                mat[i][i]=vec[i];
            }
            else{
                mat[i][j]=0;
            }
        }
    }
    return mat;
}
/*
given a matrix, return it's transpose
*/
double** transpose(double **mat, int n, int m){
    double *a;
    double **mat_t;
    int i;
    int j;
    a=calloc(n*m,sizeof(double));
    assert(a!=NULL);
    mat_t=calloc(m,sizeof(double*));
    assert(mat_t!=NULL);
    for (i=0;i<m;i++){
        mat_t[i]= a + i*n;
    }
    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            mat_t[i][j]=mat[j][i];
        }
    }
    return mat_t;
}

double** wrap(double **mat1,double **mat2, int n){
    double** mat;
    mat=sq_matrix_mul(sq_matrix_mul(transpose(mat1,n,n),mat2,n),mat1,n);
    return mat;
}
/*
given a matrix, return it's diagonal
*/
double* get_diag(double** mat, int n){
    double* vec;
    int i;
    vec=calloc(n,sizeof(double));
    assert(vec!=NULL);
    for (i=0; i<n; i++){
        vec[i]=mat[i][i];
    }
    return vec;
}

/*
given a matrix return the the cordinates of the biggest element in the matrix that isn't in the diagonal
*/
int* get_max_indexes(double** mat, int n){
    int i,j,x,y;
    double max;
    int* cor;
    cor=calloc(2,sizeof(int));
    assert(cor!=NULL);
    max=mat[0][1];
    x=0;
    y=1;
    for (i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(i!=j){
                if(abs(mat[i][j])>max){
                    max=abs(mat[i][j]);
                    x=i;
                    y=j;
                }
            }
        }
    }
    cor[0]=x;
    cor[1]=y;
    return cor;
}

double compute_theta(double x, double y, double z){
    return (x-y)/(2*z);
}

int sign(double x){
    if (x>=0){
        return 1;
    }
    else{
        return -1;
    }

}

double compute_t(double x){
    int si;
    double ab;
    double sq;

    si=sign(x);
    ab=x*si;
    sq=sqrt(pow(x,2)+1);
    return (si)/(ab+sq);
}

double compute_c(double x){
    return 1/sqrt(pow(x,2)+1);
}

double square_sum(double** mat, int n){
    double sum;
    int i,j;
    sum=0;
    for (i=0;i<n;i++){
        for (j=0;j<n;j++){
            sum+=pow(mat[i][j],2);
        }
    }
    return sum;
}

double compute_off(double** mat, int n){
    double x;
    double* diag;
    diag=get_diag(mat,n);
    x=pow(norm(diag,n),2);
    return square_sum(mat,n) - x;
}

double** many_mul(double*** mat, int n, int m){
    double **a;
    int i;
    a=mat[0];
    for (i=1;i<m;i++){
        a=sq_matrix_mul(a,mat[i],n);
    }
    return a;
}

double** get_eye_mat(int dim){
    double* arr = calloc(dim, sizeof(double));  // yoni check calloc + free
    for(int i=0;i<dim;i++){
        arr[i] = 1;
    }
    return get_diag_mat(arr, dim);
}

double** create_p_matrix(double s, double c, int x, int y, int n){
    double *a;
    double **mat;
    int i;
    int j;
    mat = get_eye_mat(n);
    mat[x][y] = -s;
    mat[y][x] = s;
    mat[x][x] = c;
    mat[y][y] = c;
    return mat;
}

double* get_deltas(double* lst, int n){
    double* n_list;
    int i;
    double x;
    n_list=calloc(n-1,sizeof(double));
    assert(n_list!=NULL);
    for (i=0; i<n-1; i++){
        x=lst[i]-lst[i+1];
        n_list[i]=x*sign(x);
    }
    return n_list;
}

int determine_k(double* lst, int n){
    int x,i,y;
    int max;
    x=0;
    y=(n+1)/2;
    max=lst[0];
    for (i=1; i<y;i++){
        if (lst[i]>max){
            max=lst[i];
            x=i;
        }
    }
    return i;
}

/*
function to swap elements
*/ 
void swap(int *a, int *b) {
  int t = *a;
  *a = *b;
  *b = t;
}

/* function to find the partition position
*/
int partition(int *array, int low, int high) {
  
  int pivot = array[high];
  int i = (low - 1);
  for (int j = low; j < high; j++) {
    if (array[j] <= pivot) {
      i++;
      swap(&array[i], &array[j]);
    }
  }
  swap(&array[i + 1], &array[high]);
  return (i + 1);
}

void quickSort(int *array, int low, int high) {
  if (low < high) {
    int pi = partition(array, low, high);
    quickSort(array, low, pi - 1);
    quickSort(array, pi + 1, high);
  }
}

int* read_file_dimensions(char* filename){
    int* arr;
    double value;
    char c;
    int vec_len;
    int vec_num;
    FILE* f = fopen(filename, "r");

    arr = calloc(2, sizeof(int));

    while (!feof(f)) {
        if(fscanf(f,"%lf%c", &value, &c) == 2){
            if(c == '\n') vec_num++;
            if(c ==',') vec_len++;
        }
     }
     fclose(f);
     arr[0] = vec_len;
     arr[1] = vec_num;
     return arr;
}


double** init_2d_array(int n, int m)
{
   //  yoni please modify this function
}

double** get_points_from_file(char* filename, int vec_len, int vec_num){
    double** points;
    FILE* f;
    int value,i,j;
    char c;


    points = init_2d_array(vec_len, vec_num);

    i=0;
    j=0;
    f = fopen(filename, "r");
    while (!feof(f)) {
        if(fscanf(f,"%lf%c", &value, &c) == 2){
            points[i][j] = value;
            if(c == '\n') i++;
            if(c ==',') j++;
        }
     }
     fclose(f);

     return points;
}


double** get_normalized_matrix(double** weighted, int dim){
    double* diag = calloc(dim, sizeof(double));  // yoni check calloc + free
    double** normalized;
    for(int i=0; i<dim; i++){
        diag[i] = sum_line(weighted, dim, i);
    }

    mul_lines(weighted, diag, dim);  // these 2 can be done in one line if we need to speed up
    mul_columns(weighted, diag, dim);

    normalized = matrix_sub(get_eye_mat(dim), weighted, dim);

    return normalized;


}

typedef struct eigen_ret{
    double** eigen_vectors;
    int k;
} EIGEN;

typedef EIGEN* EIGEN_LINK;

EIGEN_LINK get_eigens_and_k(double** normalized, int dim, int k){   // yoni please fix the memory release on this func
    int i,j;
    int* indexes;
    double temp;
    double** V;
    double** P;
    double c, s, t, theta;
    double* deltas;
    double* eigen_vals;
    EIGEN_LINK ret;
    
    int cnt =1 ;
    ret = malloc(sizeof(EIGEN));
        
    V = get_eye_mat(dim);

    print_mat(normalized, dim, dim);

            printf(" %lf \n",compute_off(normalized, dim));

    do{
        printf("\niteration nmber %d\n",cnt++);
        temp = compute_off(normalized, dim);
        indexes = get_max_indexes(normalized, dim);
        i = indexes[0];
        j = indexes[1];



        theta = compute_theta(normalized[i][i], normalized[j][j], normalized[i][j]);
        t = compute_t(theta);
        c = compute_c(t);
        s = c*t;

        printf("found vals %lf %lf %lf %lf\n", theta,t,c,s);

        P = create_p_matrix(s,c,i,j,dim);

        print_mat(P, dim, dim);
        
        printf("\n");

        V = sq_matrix_mul(V, P, dim);

        normalized = sq_matrix_mul(sq_matrix_mul(transpose(P, dim, dim) , normalized, dim ), P, dim);

        print_mat(normalized, dim, dim);

        printf(" %lf-%lf = %lf \n",temp,compute_off(normalized, dim), temp-compute_off(normalized, dim));

    } while((temp-compute_off(normalized, dim)) >= 0.000001);

    if(!k){
        eigen_vals = get_diag(normalized, dim);
        deltas = get_deltas(eigen_vals, dim);
        k = determine_k(deltas,dim-1);
    }

    ret->k = k;
    ret->eigen_vectors = V;

    printf("\nfinished jacobi\n");

    return ret;
}


int main(int argv, char* args){
    //int k=4;
    char* filename;
    char* goal;
    double** points;
    double** weighted;
    double** diagonal;
    double** normalized;
    double** eigen_vec;
    double** u_mat;
    int* dims;
    int vec_len, vec_num;
    EIGEN_LINK eigens;


/*
    dims = read_file_dimensions(filename);
    vec_len = dims[0];
    vec_num = dims[1];


    points = get_points_from_file(filename, vec_len, vec_num);

    weighted = weighted_matrix(points, vec_len, vec_num);

    normalized = get_normalized_matrix(weighted, vec_num);
    */

    double** arr = malloc(4*sizeof(double*));
    arr[0] = malloc(4*sizeof(double));
    arr[1] = malloc(4*sizeof(double));
    arr[2] = malloc(4*sizeof(double));
    arr[3]= malloc(4*sizeof(double));
    arr[0][0] = 4;
    arr[0][1]=-30;
    arr[0][2]=60;
    arr[0][3] = -35;
    arr[1][0]=-30;
    arr[1][1]=300;
    arr[1][2]=-675;
    arr[1][3] = 420;
    arr[2][0]=60;
    arr[2][1]=-675;
    arr[2][2]=1620;
    arr[2][3] = -1050;
    arr[3][0]=-35;
    arr[3][1]=420;
    arr[3][2]=-1050;
    arr[3][3] = 700;

        double** arr1 = malloc(2*sizeof(double*));
    arr1[0] = malloc(2*sizeof(double));
    arr1[1] = malloc(2*sizeof(double));

    arr1[0][0] =2;
    arr1[0][1] = 1;
    arr1[1][0] = 1;
    arr1[1][1] =2;

            double** arr2 = malloc(3*sizeof(double*));
    arr1[0] = malloc(3*sizeof(double));
    arr1[1] = malloc(3*sizeof(double));
    arr1[2] = malloc(3*sizeof(double));
    arr1[0][0] =2;
    arr1[0][1] = 1;
    arr1[0][2] = 1;
    arr1[1][0] = 1;
    arr1[1][1] =2;
    arr1[1][2] = 1;
    arr1[2][0] = 1;
    arr1[2][1] =1;
    arr1[2][2] = 2;

    vec_num = 4;
    int k=4;
    printf("gi\n");

    eigens = get_eigens_and_k(arr, vec_num, k);

    k = eigens->k;

    char sign;


    normalize_mat(eigens->eigen_vectors, vec_num, k);

    print_mat(eigens->eigen_vectors, vec_num, k);
    
    }
