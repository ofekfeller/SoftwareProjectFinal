#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

static int comp = 0;
static int mov = 0;
static int count_dis = 0;
static int hashtable[10000];

//*****************************************************************************
/* Function to sort an array using insertion sort*/
void insertionSort(int arr[], int n)
{
    int i, key, j;
    comp = 0;
    mov = 0;
    count_dis = 0;
    //int comp=0, mov=0;//comp= numbers of Comparisons, mov=nembers of moving elements
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;
        mov = mov + 1;
        /* Move elements of arr[0..i-1], that are
          greater than key, to one position ahead
          of their current position */
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
            mov = mov + 1;
            comp = comp + 2;
        }
        arr[j + 1] = key;
        mov = mov + 1;
    }

    //count distinct element
    for (int i = 0; i < n; i++)
    {
        // Move the index ahead while
        // there are duplicates
        while (i < n - 1 && arr[i] == arr[i + 1])
        {
            comp = comp + 2;
            i++;
        }
        count_dis = count_dis + 1;;
    }
    printf("\n\tInsertion Sort:");
    /*for(int i=0;i<n;i++)
    {
      printf("\t%d",arr[i]);
    }*/
    printf("\n\tnumber of distinct element: %d", count_dis);
    printf("\n\tComparisons:%d\tMoving:%d", comp, mov);
}

//*****************************************************************************
// Counting sort in C programming
int  counting_sort(int a[], int n)
{
    comp = 0;
    mov = 0;
    count_dis = 0;
    int i, j;
    int b[101], c[101];
    for (i = 0; i <= 100; i++) {
        c[i] = 0;
        mov = mov + 1;
    }
    printf("a");
    for (j = 0; j < n; j++) {
        c[a[j]] = c[a[j]] + 1;
        mov = mov + 1;
    }
    printf("b");
    //count distinct element
    for (int i = 0; i <= 100; i++)
    {
        comp = comp + 1;
        if (c[i] != 0) {
            count_dis = count_dis + 1;
        }
    }

    printf("c");
    for (i = 1; i < 101; i++) {
        c[i] = c[i] + c[i - 1];
        mov = mov + 1;
    }
    printf("d");
    for (j = n - 1; j >= 0; j--)
    {
        printf("%d", j);
        b[c[a[j]]] = a[j];
        c[a[j]] = c[a[j]] - 1;
        mov = mov + 2;
    }


    printf("\n\tCounting Sort:");
    /*for(int i=1;i<=n;i++)
    {
      printf("\t%d",b[i]);
    }*/
    printf("\n\tnumber of distinct element: %d", count_dis);
    printf("\n\tComparisons:%d\tMoving:%d", comp, mov);
}


// Heap Sort in C*****************************************************

  // Function to swap the the position of two elements
void swap(int* a, int* b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void heapify(int arr[], int n, int i) {
    // Find largest among root, left child and right child
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;

    // If left child is larger than root
    if (left < n && arr[left] > arr[largest]) {
        largest = left;
        comp = comp + 1;
    }

    // If right child is larger than largest so far
    if (right < n && arr[right] > arr[largest]) {
        largest = right;
        comp = comp + 1;
    }
    // Swap and continue heapifying if root is not largest
    if (largest != i) {
        mov = mov + 1;
        swap(&arr[i], &arr[largest]);
        heapify(arr, n, largest);
    }
}

// Main function to do heap sort
void heapSort(int arr[], int n) {
    comp = 0;
    mov = 0;
    count_dis = 0;
    //int comp=0, mov=0;//comp= numbers of Comparisons, mov=nembers of moving elements

    // Build max heap
    for (int i = n / 2 - 1; i >= 0; i--){
        heapify(arr, n, i);
    }
    printf("a");
    // Heap sort
    for (int i = n - 1; i >= 0; i--) {
        swap(&arr[0], &arr[i]);
        mov = mov + 1;
        // Heapify root element to get highest element at root again
        heapify(arr, i, 0);
    }
    printf("b");
    //count distinct element
    for (int i = 0; i < n; i++)
    {
        // Move the index ahead while
        // there are duplicates
        while (i < n - 1 && arr[i] == arr[i + 1])
        {
            i++;
            comp = comp + 1;
        }
        count_dis = count_dis + 1;;
    }
    printf("c");
    printf("\n\tHeap Sort:");
    /*for(int i=0;i<n;i++)
    {
      printf("\t%d",arr[i]);
    }*/
    printf("\n\tnumber of distinct element: %d", count_dis);
    printf("\n\tComparisons:%d\tMoving:%d", comp, mov);
}




//***************************************************

int search(int value)
{
    comp = comp + 1;
    if (hashtable[value] == value)
        return 1;
    else
        return 0;
}
void f_hashtable(int arr[], int n)
{
    comp = 0;
    mov = 0;
    count_dis = 0;
    // Creates an empty hashset
    for (int i = 0; i < n; i++) {
        printf("%d\n",i);
        hashtable[i] = -1;
        mov = mov + 1;
    }
    printf("a");
    // Traverse the input array
    for (int i = 0; i < n; i++) {

        // If not present, then put it in
        // hashtable and increment result
        //0=the key not present
        //1= the key present
        if (search(arr[i]) == 0) {
            hashtable[arr[i]] = arr[i];
            mov = mov + 1;
            count_dis = count_dis + 1;

        }
        comp = comp + 1;
    }
    printf("b");
    printf("\n\tHash table:\n\tnumber of distinct element: %d", count_dis);
    printf("\n\tComparisons:%d\tMoving:%d", comp, mov);
}

//***********************************************************************
void D(int A[], int n) { // A is an array of numbers
    comp = 0;
    mov = 0;
    int U_Size = 1;
    bool U;
    for (int i = 1; i < n; i++)
    {
        U = true;
        for (int j = 0; j <= U_Size; j++)
        {
            comp = comp + 1;
            if (A[j] == A[i])
            {
                U = false;
                j = U_Size + 1;
            }
        }
        comp = comp + 1;
        if (U == true)
        {
            U_Size = U_Size + 1;
            A[U_Size] = A[i];
            mov = mov + 1;
        }
    }
    printf("\n\tWith source algoritem");
    /*for(int i=0;i<n;i++)
    {
      printf("\t%d",A[i]);
    }*/
    printf("\n\tnumber of distinct element: %d", U_Size + 1);
    printf("\n\tComparisons:%d\tMoving:%d", comp, mov);
}



//duplicate array
void duplicateArray(int source[], int dest[], int size)
{
    for (int i = 0; i < size; i++)
    {
        dest[i] = source[i];
    }
}

int main(void)
{
    int size;
    printf("pls enter array size \n");
    scanf("%d", &size);

    int* arr = malloc(sizeof(int) * size);
    int* arr1 = malloc(sizeof(int) * size);
    int* arr2 = malloc(sizeof(int) * size);
    int* arr3 = malloc(sizeof(int) * size);
    int* arr4 = malloc(sizeof(int) * size);
    int* arr5 = malloc(sizeof(int) * size);

    for (int i = 0; i < size; i++) {
        arr[i] = rand() % 100 + 1;
        //printf("%d\n", arr[i]);
    }
    
        
    
      /*printf("\nElements of the array:");
      for(i=0;i<sz;i++)
      {
        printf("\t%d",arr[i]);
      }*/

    duplicateArray(arr, arr1, size);
    duplicateArray(arr, arr2, size);
    duplicateArray(arr, arr3, size);
    duplicateArray(arr, arr4, size);
    duplicateArray(arr, arr5, size);

    printf("\n\nfor N=%d:", size);
    f_hashtable(arr5, size);//מימוש סעיף 5 - טבלת גיבוב
    D(arr1, size); // מימוש סעיף 1 - אלגוריתם מקורי
    free(arr1);
    insertionSort(arr2, size);//מימוש סעיף 2 - מיון הכנסה
    free(arr2);
    heapSort(arr3, size);// מימוש סעיף 3 - מיון מבוסס השוואות 
    free(arr3);


    
    free(arr);



    free(arr4);
    free(arr5);

    return 0;
}