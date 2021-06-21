#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double norm(double *a, double *b, int size);
int classify(double *vec, double **centroids , int size_vec, int num_cent);
double * update_center(INDEX_LINK a, double **points, int vect_size);

struct index_linked_list{
    int data;
    struct index_linked_list *next
};

typedef struct index_linked_list INDEX;
typedef INDEX* INDEX_LINK;

//fix, first element doesnt contain data
double * update_center(INDEX_LINK a, double **points, int vect_size){
    double *vec;
    vec=(double*)malloc(vect_size * sizeof(double));
    int count=0;
    int c=0;
    while(a != NULL){
        count++;
        c=a->data;
        for (int j=0;j<vect_size;j++){
            vec[j]+=points[c][j];
        }
        a=a->next;
    }
    for (int i=0;i<vect_size;i++){
        vec[i]=vec[i]/count;
    }
    return vec;
}

double norm(double *a,double *b, int size){
    double norm=0;
    double c=0;
    double *vec;
    vec=(double*)malloc(size * sizeof(double));
    for (int i=0;i<size; i++){
        c = (b[i]-a[i]);
        vec[i]=c*c;
    }
    for (int i=0; i<size; i++){
        norm+=vec[i];
    }
    return norm;
}

int classify(double *vec, double **centroids, int size_vec, int num_cent){
    int min_ind=0;
    double norm1=0;
    int min_norm=norm(vec,centroids[0],size_vec);
    for (int i=1;i<num_cent;i++){
        norm1=norm(vec,centroids[i],size_vec);
        if (norm1<min_norm){
            min_norm=norm1;
            min_ind=i;
        }
    }
    return min_ind;
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

void remove_from_old(int i, INDEX_LINK head){
    int c=0;
    INDEX_LINK temp=NULL;
    temp = (INDEX*)malloc(sizeof(INDEX));
    c=head->data;
    if (i==c){
        head->next=head->next->next;
    }
    else{
        temp=head->next;
        while (temp!=NULL){
            c=temp->data;
            if(c==i){
                head->next=temp->next;
            }
            temp=temp->next;
            head=head->next;
        }
    }
}