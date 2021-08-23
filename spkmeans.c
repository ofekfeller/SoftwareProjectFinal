#include <Math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

double epsilon = 0.001;

int classify(double *vec,int num_cent, int size_vec, double centers[][size_vec]){
    int min_ind=0;
    double tmp_norm=0;
    double min_norm;
    int i;

    min_norm = norm(sub_vectors(vec,centers[0],size_vec),size_vec);


    for (i=1; i<num_cent; i++){

        tmp_norm = norm(sub_vectors(vec,centers[i],size_vec),size_vec);
        tmp_norm = sec_norm(vec,centers[i],size_vec); 

        if (tmp_norm<min_norm){
            min_norm=tmp_norm;
            min_ind=i;
        }
    }
    return min_ind;
}

void copy_array(double *source, double *new, int dim){
    int i;
    for (i = 0; i < dim; i++)
    {
        source[i] = new[i];
    }
}

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
    }

    return clusters;
}

void update_center(double* origin, CEN_INFO cluster, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        origin[i] = (cluster.sum[i] / cluster.cnt);
    }
    
}

double* sub_vectors(double *vec1,double *vec2, int n){
    int i;
    double* vec=calloc(n,sizeof(double));
    assert(vec!=NULL);
    for (i=0; i<n;i++){
        vec[i]=vec1[i]-vec2[i];
    }
    return vec;
}

void add_vectors(double *vec1, double *vec2, int n){
    int i;
    double tmp;

    for (i = 0; i < n; i++)
    {
        tmp = vec1[i] + vec2[i];
        vec1[i] = tmp;  
    }
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

double get_abs(double num){
    if(num<0){
        return -1*num;
    }
    return num;
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
                if(get_abs(mat[i][j])>max){
                    max=get_abs(mat[i][j]);
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
    double ret;
    ret = x-y;
    ret = ret / (2*z);
    return ret;
}

double sign(double x){
    if (x>=-0.000000000001){
        return 1;
    }
    else{
        return -1;
    }

}

double compute_t(double x){
    double si;
    double sq;

    si=sign(x);
    sq=sqrt(pow(x,2)+1);
    return (si)/(get_abs(x)+sq);
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
    double** mat;
    double* arr = calloc(dim, sizeof(double));  // yoni check calloc + free
    for(int i=0;i<dim;i++){
        arr[i] = 1;
    }
    //split the return
    //return get_diag_mat(arr, dim);
    mat= get_diag_mat(arr, dim);
    free(arr);
    return mat;
}

double** create_p_matrix(double s, double c, int x, int y, int n){
    double *a;
    double **mat;
    int i;
    int j;
    mat = get_eye_mat(n);
    mat[x][y] = s;
    mat[y][x] = -s;
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
   double *a;
    double **mat;
    int i;
    int j;
    a=calloc(n*m,sizeof(double));
    assert(a!=NULL);
    mat=calloc(n,sizeof(double*));
    assert(mat!=NULL);
    for (i=0;i<n;i++){
        mat[i]= a + i*m;
    }
    return mat;
}

double** get_points_from_file(char* filename, int vec_len, int vec_num){
    double** points;
    FILE* f;
    int value,i,j;
    char c;


    points = init_2d_array(vec_num, vec_len);

    i=0;
    j=0;
    f = fopen(filename, "r");
    while (!feof(f)) {
        if(fscanf(f,"%lf%c", &value, &c) == 2){
            points[i][j] = value;
            if(c == '\n'){
                i++;
                j=0;
            }
            if(c ==',') j++;
        }
     }
     fclose(f);

     return points;
}

double* get_diag_vec(double** weighted, int dim){
    double* diag = calloc(dim, sizeof(double));  // yoni check calloc + free
    
    for(int i=0; i<dim; i++){
        diag[i] = sum_line(weighted, dim, i);
    }
    return diag;
}


double** get_normalized_matrix(double** weighted, double* diag, int dim){ // yoni check calloc + free
    double** normalized;

    mul_lines(weighted, diag, dim);  // these 2 can be done in one line if we need to speed up
    mul_columns(weighted, diag, dim);

    normalized = matrix_sub(get_eye_mat(dim), weighted, dim);
    free(diag);
    return normalized;


}

double** deep_copy(double** mat, int dim){
    double** ret;

    ret = init_2d_array(dim,dim);

    for(int i = 0;i<dim;i++){
        for(int j=0;j<dim;j++){
            ret[i][j] = mat[i][j];
        }
    }
    return ret;

}

void compute_normalized(double** mat, int dim, double c, double s, int i, int j){

    // handle memory issues


    double temp;
    double** copy;

    copy = deep_copy(mat,dim);
    for(int k=0;k<dim;k++){
        if(k == i || k == j){
        continue;}
        temp = c*copy[k][i]-s*copy[k][j];
        mat[k][i] = temp;
        mat[i][k] = temp;
    }

    for(int t=0;t<dim;t++){
        if(t == i || t == j){
        continue;}
        temp = c*copy[t][j]+s*copy[t][i];
        mat[t][j] = temp;
        mat[j][t] = temp;
    }


    mat[i][i] = (pow(c,2)*copy[i][i])+(pow(s,2)*copy[j][j])-(2*s*c*copy[i][j]);  
    mat[j][j] = (pow(c,2)*copy[j][j])+(pow(s,2)*copy[i][i])+(2*s*c*copy[i][j]); 

    mat[i][j] = 0;
    mat[j][i] = 0;

}

typedef struct eigen_ret{
    double * eigen_values;
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
    double** t_P;
    double c, s, t, theta;
    double* deltas;
    double* eigen_vals;
    EIGEN_LINK ret;
    
    int cnt =1 ;
    ret = malloc(sizeof(EIGEN));
        
    V = get_eye_mat(dim);

    do{
        temp = compute_off(normalized, dim);
        indexes = get_max_indexes(normalized, dim);
        i = indexes[0];
        j = indexes[1];

        theta = compute_theta(normalized[j][j], normalized[i][i], normalized[i][j]);
        t = compute_t(theta);
        c = compute_c(t);
        s = t*c;

        P = create_p_matrix(s,c,i,j,dim);

        V = sq_matrix_mul(V, P, dim);    // check what is the formula to compute P

        compute_normalized(normalized,dim,c,s,i,j);

    } while((temp-compute_off(normalized, dim)) > epsilon);

    if(!k){
        eigen_vals = get_diag(normalized, dim);
        deltas = get_deltas(eigen_vals, dim);
        k = determine_k(deltas,dim-1);
    }

    ret->k = k;
    ret->eigen_vectors = V;
    ret->eigen_values = get_diag(normalized, dim);
    
    free(indexes);
    free(V);
    free(P);
    free(eigen_vals);
    free(deltas);
    return ret;
}

void kmeans_goal(double** points, char* goal, int vec_num, int dim){
    EIGEN_LINK eigens;
    double** weighted;
    double** normalized;
    double* diag;

    if(goal=="jacobi"){
        eigens = get_eigens_and_k(points, vec_num, vec_num);
        print_vec(eigens->eigen_values, vec_num);
        print_mat(transpose(eigens->eigen_vectors, vec_num, vec_num), vec_num, vec_num);

    }

    weighted = weighted_matrix(points, dim, vec_num);

    if(goal=="wam"){
        print_mat(weighted, vec_num, vec_num);
    }

    diag = get_diag_vec(weighted, vec_num);

    if(goal=="ddg"){
        print_mat(get_diag_mat(diag, vec_num), vec_num, vec_num);
    }

    normalized = get_normalized_matrix(weighted,diag, vec_num);

    if(goal=="lnorm"){
        print_mat(normalized, vec_num, vec_num);
    }

    free(normalized);
    free(diag);
    free(weighted);
    free(eigens);
}

double** get_spk_points(double** points, int dim, int vec_num, int k){
    EIGEN_LINK eigens;
    double** weighted;
    double** normalized;
    double* diag;
    double** centers;
    CEN_LINK clusters;
    int max_iter=300;
    weighted = weighted_matrix(points, dim, vec_num);
    diag = get_diag_vec(weighted, vec_num);
    normalized = get_normalized_matrix(weighted,diag, vec_num);
    eigens = get_eigens_and_k(normalized, vec_num, k);
    centers=init_2d_array(k,k);
        
        for(i=0;i<k;i++){
            copy_array(centers[i], eigens->eigen_vectors[i], k);
        }
        int vec_to_cen[vec_num];
        CEN_LINK clusters = init_clusters(k, k);

        int center;

        for (i = 0; i < vec_num; i++)
        {
        center = classify(eigens->eigen_vectors[i],  k, k, centers);

        vec_to_cen[i] = center;

        add_vectors(clusters[center].sum , eigens->eigen_vectors[i], k);

        clusters[center].cnt++;

        }

        int bool;
        int old_center;

        bool = 1;
        for (int i=0; i<max_iter-1;i++){

            for(int j=0;j<vec_num;j++){

                center = classify(eigens->eigen_vectors[j] ,k, k, centers);

                old_center = vec_to_cen[j];

                if (old_center!=center){

                    sub_vectors(clusters[old_center].sum, eigens->eigen_vectors[j], k);

                    add_vectors(clusters[center].sum , eigens->eigen_vectors[j], k);

                    clusters[center].cnt++;
                    clusters[old_center].cnt--;

                    vec_to_cen[j] = center;
                    
                    bool = 1;
                }

            }
        
            if (bool == 0){
                
                break;
                
            }  
            
            for (int j = 0; j < k; j++)
            {
                
                update_center(centers[j], clusters[j], k);
                
            }

            bool = 0;
        }
    free(normalized);
    free(diag);
    free(weighted);
    free(eigens);
    free(clusters);
    free(vec_to_cen);
    return centers;
}

int main(int argv, char* args){
    int k;
    double** points;
    double** weighted;
    double** diagonal;
    double** normalized;
    double** eigen_vec;
    double** u_mat;
    int* dims;
    int vec_len, vec_num;
    EIGEN_LINK eigens;
    char* goal;
    char* file_name;
    double* diag;

   assert(argv==4);
    k = (int)args[0];
    goal = args[1];
    file_name = args[2];

    dims = read_file_dimensions(file_name);
    vec_len = dims[0];
    vec_num = dims[1];


    points = get_points_from_file(file_name, vec_num, vec_len);

    if(goal=="jacobi"){
        eigens = get_eigens_and_k(points, vec_num, vec_num);
        print_vec(eigens->eigen_values, vec_num);
        print_mat(transpose(eigens->eigen_vectors, vec_num, vec_num), vec_num, vec_num);
        return 1;

    }

    weighted = weighted_matrix(points, vec_len, vec_num);

    if(goal=="wam"){
        print_mat(weighted, vec_num, vec_num);
        return 1;
    }

    diag = get_diag_vec(weighted, vec_num);

    if(goal=="ddg"){
        print_mat(get_diag_mat(diag, vec_num), vec_num, vec_num);
        return 1;
    }

    normalized = get_normalized_matrix(weighted,diag, vec_num);

    if(goal=="lnorm"){
        print_mat(normalized, vec_num, vec_num);
        return 1;
    }

    /*
    double** arr = malloc(4*sizeof(double*));
    arr[0] = malloc(4*sizeof(double));
    arr[1] = malloc(4*sizeof(double));
    arr[2] = malloc(4*sizeof(double));
    arr[3]= malloc(4*sizeof(double));
    arr[0][0] = 4;
    arr[0][1] = -30;
    arr[0][2] = 60;
    arr[0][3] = -35;
    arr[1][0] = -30;
    arr[1][1] = 300;
    arr[1][2] = -675;
    arr[1][3] = 420;
    arr[2][0] = 60;
    arr[2][1] = -675;
    arr[2][2] = 1620;
    arr[2][3] = -1050;
    arr[3][0] = -35;
    arr[3][1] = 420;
    arr[3][2] = -1050;
    arr[3][3] = 700;

    double** arr1 = malloc(2*sizeof(double*));
    arr1[0] = malloc(2*sizeof(double));
    arr1[1] = malloc(2*sizeof(double));

    arr1[0][0] =2;
    arr1[0][1] = 1;
    arr1[1][0] = 1;
    arr1[1][1] =2;

    double** arr2 = malloc(3*sizeof(double*));
    arr2[0] = malloc(3*sizeof(double));
    arr2[1] = malloc(3*sizeof(double));
    arr2[2] = malloc(3*sizeof(double));
    arr2[0][0] =1;
    arr2[0][1] = sqrt(2);
    arr2[0][2] = 2;
    arr2[1][0] = sqrt(2);
    arr2[1][1] =3;
    arr2[1][2] = sqrt(2);
    arr2[2][0] = 2;
    arr2[2][1] =sqrt(2);
    arr2[2][2] = 1;

    vec_num = 4;
    k=4;
    */

    eigens = get_eigens_and_k(normalized, vec_num, 0);
    normalize_mat(eigens->eigen_vectors, vec_num, k);

    k = eigens->k;   

    printf("%d\n",k); 

    print_vec(eigens->eigen_values, vec_num);

    print_mat(eigens->eigen_vectors, vec_num, k);

    k_means(eigens->eigen_vectors,vec_num, k);

    }
