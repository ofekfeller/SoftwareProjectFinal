#include <stdio.h>
#include <stdlib.h>


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

    clusters = calloc(k, sizeof(CEN_INFO));
    for ( i = 1; i < k; i++)
    {
        clusters[i].cnt=0;
        clusters[i].sum = calloc(dim, sizeof(double));
    }

    return clusters;
}


double* main(){

    int k = 3; //change this
    int max_iter= 200; //change this 
    int i;
    int j;

    double value;
    char c;
    ELEMENT vec;
    LINK tail;
    LINK head;
    int flag = 1;
    int counter = 0;
    if (scanf("%lf%c", &value, &c) == 2){
        head = (ELEMENT*)malloc(sizeof(ELEMENT));
        head->data=value;
        tail=head;
        counter++;
    }
    while (scanf("%lf%c", &value, &c)==2){
        tail->next = (ELEMENT*)malloc(sizeof(ELEMENT));
        tail = tail->next;
        tail->data=value;
        counter++;
        if (c=='\n'){
            tail->next = NULL;
            break;
        }
    }

    double *arr = calloc(counter,sizeof(double));
    for (i=0; i<counter;i++){
        arr[i]=head->data;
        head=head->next;
    }
    MAIN_LINK points_head = NULL,points_tail=NULL;
    points_head = (MAIN_ELEMENT*)malloc(sizeof(MAIN_ELEMENT));
    points_head->vec=arr;
    points_tail = points_head;
    double *vector = calloc(counter,sizeof(double));
    flag = 0;
    int vec_cnt = 0;
    while (scanf("%lf%c", &value,&c) == 2)
    {
        vector[flag] = value;
        flag++;
        if (c == '\n')
        {
        points_tail->next = (MAIN_ELEMENT*)malloc(sizeof(MAIN_ELEMENT));
        points_tail = points_tail->next;
        points_tail->vec=vector;
        flag = 0;
        vector = (double*)malloc(counter*sizeof(double)); 
        vec_cnt++;
        } 

    }
    tail->next= NULL;

    // creating the final 2-d array of points

    double points[vec_cnt][counter];
    MAIN_LINK pointer = points_head;
    for(int i=0;i<vec_cnt;i++){
        
        copy_array(points[i], pointer->vec, counter); 

        pointer = pointer->next;

    }

    // initializing k centers
    double centers[k][counter];

    for(i=0;i<k;i++){

        copy_array(centers[i], points[i], counter);

    }

    char sign;
    for (i=0; i<k; i++){

       print_array(centers[i],counter);
       printf("\n");

    }

    printf("Creating clusters\n");

    int vec_to_cen[vec_cnt];

    CEN_LINK clusters = init_clusters(k, vec_cnt);
    int center;

    for (i = 0; i < vec_cnt; i++)
    {
       printf("%d\n", i);

       center = classify(points[i], centers, counter, k);
       vec_to_cen[i] = center;
       add_array(clusters[center].sum , points[i], counter);
       clusters[center].cnt++;

    }

    printf("starint iters clusters");

    int bool;
    int old_center;
    for (int i=0; i<max_iter-1;i++){

        printf("finished %d", i);

        bool = 0;

        for(int j=0;j<vec_cnt;j++){

            center = classify(points[j] ,centers, counter, k);

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

        for (int i = 0; i < k; i++)
        {
            update_center(centers[k], clusters[k], counter);
        }
    }

    printf("finished Creating clusters");

    for (int i=0; i<k; i++){

        for (int j; j < counter; j++)
        {
            if (j==counter-1) sign='\n';
            else sign = ',';
            printf("%lf%c", centers[i][j], sign);
        }
    }
}

void copy_array(double *source, double *new, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        source[i] = new[i];
    }
    
}

void add_array(double *source, double *new, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        source[i] += new[i];
    }
    
}

void sub_array(double *source, double *new, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        source[i] -= new[i];
    }
    
}

void update_center(double* origin, CEN_INFO cluster, int dim){
    int i;

    for (i = 0; i < dim; i++)
    {
        origin[i] = (cluster.sum[i] / dim);
    }
    
}


int classify(double *vec,int num_cent, int size_vec, double **centers){
    int min_ind=0;
    double tmp_norm=0;
    double min_norm;

    printf("try %lf \n", centers[0][1]);

    print_array(vec,size_vec);

    print_array(centers[0],size_vec);

    min_norm = norm(vec, centers[0], size_vec);


    for (int i=1;i<num_cent;i++){
        tmp_norm = norm(vec,centers[i],size_vec); 

        if (tmp_norm<min_norm){
            min_norm=tmp_norm;
            min_ind=i;
        }
    }
    return min_ind;
}

double norm(double *array1,double *array2, int size_1){
    double sum1;
    double tmp_norm;
    int i;

    sum1=0;
    tmp_norm=0;

    printf("start norm \n");

    for (i = 0; i < size_1; i++)
    {
        printf("hi");
    }
    

    printf("finish norm");
    return tmp_norm;
}

int update_row(double * to, double * from, int k){
    // 1 - changed, 0 - false
    int cnt=0;
    for (int i=0; i<k;i++){
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

void print_array(double *array, int len){
    int i;

    for ( i = 0; i < len; i++)
    {
        printf("%lf,",array[i]);
    }
}
