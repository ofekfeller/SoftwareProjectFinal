#define G(y) (y+y*y)

/*void main(int argc, char* argv[]){
    int flag =argc-1;
    int i=0;
    while (flag!=0){
        if (strlen(argv[flag])==1){
            i++;
        }
        flag--;
    }
    printf("%d",i);
}



void concat(LIST* node1, LIST* node2){
    while(node1.next != NULL){
        node1=node1.next;
    }
    node1.next = node2;
}


int** func(int n){
    int** res;
    res = (int*)calloc(n, sizeof(int*));
    for(int i=0; i<n; i++){
        res[i] = calloc(n, sizeof(int));
    }

    return res;
}

void main(){
    char str1[]="ofek";
    char str2[]="feller";
    printf("%d",strlen(str1));
    printf("%d",strlen(str2));
    int ner_len = strlen(str1) + strlen(str2);
    char* res = (char*)malloc(ner_len*sizeof(char));
    int i=0;

    while(str1[i] != '\0'){
        res[2*i] = str1[i];
        res[2*i+1] = str2[i];
        i++;
    }
    res[2*i]='\0';
    printf("%s", res);
}


void remove(LINK* ar1, LINK* ar2){
    LINK* pointer1 = ar1.next;
    LINK* pointer2 = ar2.next;
    LINK* prev = arr2;
    LINK* temp;
    while(pointer1 !=NULL && pointer2!=NULL){
        if(pointer1.value == pointer2.value){
            pointer2 = remove_node(prev, pointer2);
        }
        elif (pointer1.value > pointer2.value){
            prev = pointer2;
            pointer2= pointer2.next;
        }
        else{
            pointer1 = pointer1.next;
        }
    }
}

void remove_node(LINK* first, LINK* sec){
    first.next = sec.next;
    free(sec);
    return first.next;
}


void main(int argc, char const *argv[])
{
    int n=6;
    char** res = (char**)malloc(n*sizeof(char*));
    for (int i=0;i<n;i++){
        res[i] = malloc((i+1)*sizeof(char));
        for (int j = 0; j < i; j++)
        {
            res[i][j]='#';
        }
        res[i][i] = '\0';
    }

    for(int i=0; i<n;i++){
        printf("%s\n", res[i]);
    }
}


char* duplicate(char* str){
    int num = strlen(str);
    char* res = (char*)malloc((2*num+1)*sizeof(char));
    for(int i=0; i<num;i++){
        res[i]=str[i];
        res[num+i]=str[i];
    }
    res[2*num] = '\0';
    return res;
}

void main()
{
    /*
    assert(argc>1);
    int n = atoi(argv[1]);
    int counter = 0;
    while(n>0){
        counter++;
        n= n/10;
    }
    return counter;
 
   char ofek[]="ofek";
   
    printf("%c",'c'+30);
   //printf("%s", duplicate(ofek));
}




int eq(char* str1, char* str2){
    if(strlen(str1)!=strlen(str2)){
        return 0;
    }
    for (int i = 0; i < strlen(str1); i++)
    {
        if(str1[i]!=str2[i]) return 0;
    }
    return 1;
}

typedef struct node{
    int val;
    NODE* next;
}NODE;

typedef struct list{
    NODE* head;
    NODE* tail;
}LIST;


LIST *Bucket(void *A, int nElementsA, int sizeOfAnElementA, void *options, int nOptions, int (*comp_func)(void*, void*)){
    LIST** res = (LIST**)malloc(nOptions*sizeof(LIST*));
    for(int i=0; i<nOptions; i++){
        res[i] = (LIST*)malloc(sizeof(LIST));
        for (int j = 0; j < sizeOfAnElementA; j++)
        {
            if(comp_func(options+i, A+j))
            {
                insert_node(res[0], (&A)+j);
            }   
        }
        
    }

return res;
}


void insert_node(LIST* list1, int new_val){
    if(list1->head == NULL){
        list1->head = (NODE*)malloc(sizeof(NODE));
        list1->head->val = new_val;
        list1->head->next = NULL;
        list1->tail = list1->tail;
    }
    else{
        NODE* new_node = (NODE*)malloc(sizeof(NODE));
        new_node->val = new_val;
    
        list1->tail->next = new_node;
        list1->tail = new_node;        
    }
}
*/
void main(){
    char* s1 = "9099hel1lo7";
    printf("%d", atoi(s1));
}